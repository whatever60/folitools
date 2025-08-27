from pathlib import Path

import pandas as pd
import polars as pl
import networkx as nx

from .utils import expand_path_to_list
from .primer_selection.recover_utils import simplify_gene_list

expected = [
    "read_id",
    "contig",
    "position",
    "gene",
    "umi",
    "umi_count",
    "final_umi",
    "final_umi_count",
    "unique_id",
]


# def deduplicate_umi(
#     group_count: pl.DataFrame | pl.LazyFrame,
# ) -> pd.DataFrame:
#     gene_counts = (
#         group_count.select(["gene", "final_umi"])
#         .unique()
#         .group_by("gene")
#         .agg(pl.len().alias("umi_count"))
#         .select(["gene", "umi_count"])
#     )
#     if isinstance(gene_counts, pl.LazyFrame):
#         gene_counts = gene_counts.collect()
#     gene_counts = gene_counts.to_pandas().set_index("gene")
#     return gene_counts


def count_umi(group_count: pl.DataFrame | pl.LazyFrame) -> pd.DataFrame:
    result = (
        normalize_multimapping_umitools_group_output(group_count.lazy())
        .select(pl.struct(["gene", "final_umi"]).value_counts())
        .select(
            umi=pl.col("gene").struct.field("gene").struct.field("gene")
            + ":"
            + pl.col("gene").struct.field("gene").struct.field("final_umi"),
            count=pl.col("gene").struct.field("count"),
        )
    )
    if isinstance(result, pl.LazyFrame):
        result = result.collect()
    return result.to_pandas().set_index("umi")


def count_gene(group_count: pl.DataFrame | pl.LazyFrame) -> pd.DataFrame:
    result = (
        normalize_multimapping_umitools_group_output(group_count.lazy())
        .select(pl.col("gene").value_counts())
        .select(
            gene=pl.col("gene").struct.field("gene"),
            count=pl.col("gene").struct.field("count"),
        )
    )
    if isinstance(result, pl.LazyFrame):
        result = result.collect()
    return result.to_pandas().set_index("gene")


def consolidate_read_id(df: pl.LazyFrame, separator: str = ",") -> pl.LazyFrame:
    """Build an undirected gene-gene co-occurrence edge list from a tall table.

    Each read is treated as a "bag of genes": we first collect the set of unique
    genes per read (collapsing multiple rows of the same read), then add an
    undirected edge for every unordered pair of genes that co-occur in that read.
    The final edge weight is the number of distinct reads supporting that pair.

    Args:
        df: Input table with at least ``id_col`` and ``gene_col``.
        id_col: Column identifying a read (non-unique; will be grouped).
        gene_col: Column containing one or more gene IDs; multiple genes may be
            stored in a single cell, separated by ``separator``.
        separator: Separator used inside ``gene_col`` for multiple genes.

    Returns:
        A Polars DataFrame with columns:
            - ``gene_a`` (str): One endpoint of the edge.
            - ``gene_b`` (str): The other endpoint (lexicographically > ``gene_a``).
            - ``weight`` (i64): Number of distinct reads where the pair co-occurs.

    Notes:
        - Empty / null genes are dropped.
        - Duplicate genes within the same read are deduplicated before pairing.
        - Edges are undirected; we enforce a consistent ordering (a < b).
    """
    # 1) Normalize to one gene per row, keep only (read, gene) unique combos.
    read_gene = (
        df.with_columns(
            pl.col("gene")
            .cast(pl.Utf8)
            .str.split(separator)
            .list.eval(pl.element().str.strip_chars())
            .alias("_gene_list")
        )
        .explode("_gene_list")
        .with_columns(pl.col("_gene_list").alias("gene"))
        .filter(pl.col("gene").is_not_null() & (pl.col("gene") != ""))
        .unique(["read_id", "gene"])
    )

    # 2) Collect unique genes per read
    per_read = read_gene.group_by("read_id").agg(
        [pl.col("gene").unique().sort(), pl.col("sample").first()]
        + [pl.col(c).first() for c in expected if c != "gene" and c != "read_id"]
    )
    return per_read


def per_read_to_edges(df: pl.LazyFrame, id_col: str = "read_id") -> pl.LazyFrame:
    """Convert per-read gene lists to edge list."""
    # 3) Create all unordered pairs within each read via a self-join on id_col.
    #    Because each read's gene list is unique and sorted, filtering gene_a < gene_b
    #    yields each pair exactly once per read.
    long = df.select(id_col, pl.col("gene")).explode("gene").rename({"gene": "gene_a"})
    pairs_per_read = (
        long.join(long, on=id_col, suffix="_b")
        .filter(pl.col("gene_a") < pl.col("gene_a_b"))
        .select(
            pl.col("gene_a"),
            pl.col("gene_a_b").alias("gene_b"),
        )
    )

    # 4) Collapse across reads: count distinct reads supporting each pair.
    edges = pairs_per_read.group_by(["gene_a", "gene_b"]).len().rename({"len": "count"})

    return edges.sort("count", descending=True)


def edges_to_graph(
    edges: pl.LazyFrame,
    *,
    gene_a_col: str = "gene_a",
    gene_b_col: str = "gene_b",
    weight_col: str = "count",
    all_genes: list[str] | None = None,
) -> nx.Graph:
    """Convert an undirected, weighted edge list (Polars) to a NetworkX graph.

    Args:
        edges: Polars DataFrame with at least two endpoint columns and one weight column.
        gene_a_col: Column name for endpoint A.
        gene_b_col: Column name for endpoint B.
        weight_col: Column name for edge weight (e.g., supporting read count).
        all_genes: Optional list of all gene IDs to include as nodes. If provided,
            genes with no edges (singletons) will be present as isolated nodes.
        min_weight: Keep only edges whose weight is >= this threshold.

    Returns:
        A NetworkX undirected Graph with 'weight' on edges.
    """
    g = nx.Graph()
    if all_genes is not None:
        g.add_nodes_from(all_genes)

    # Add weighted edges
    for a, b, w in (
        edges.select(gene_a_col, gene_b_col, weight_col).collect().iter_rows()
    ):
        g.add_edge(a, b, weight=int(w))

    return g


def process_count_file(fps: list[str], dedup_umi: bool = True) -> pd.DataFrame:
    """
    Process a single count file—either raw UMI records or a precomputed gene to count table.

    Args:
        fp: Path to a tab-delimited count file.
            If the header matches
            ["read_id", "contig", "position", "gene", "umi",
             "umi_count", "final_umi", "final_umi_count", "unique_id"]
            it will compute per-gene UMI counts. Otherwise it expects
            a two-column table (gene and count).

    Returns:
        pd.Series: Index is gene, values are umi_count, name is sample (file stem).
    """
    series_list = []
    samples = []
    for fp in fps:
        sample = Path(fp).stem.split(".")[0]

        header = pd.read_csv(fp, nrows=0, sep="\t").columns.tolist()

        if header == expected:  # A table generated by `umi_tools group`
            series_list.append(
                pl.scan_csv(fp, separator="\t")
                .filter(pl.col("gene") != "Unassigned")
                .with_columns(sample=pl.lit(sample))
            )
        else:  # A table generated by `umi_tools count`
            df = pd.read_csv(fp, sep="\t", index_col=0)
            if df.shape[1] != 1:
                raise ValueError(f"Expected 1 column in {fp!r}, got {df.shape[1]}")
            s = df.iloc[:, 0]
            s.name = sample
            series_list.append(s)
        samples.append(sample)

    lazy_dfs = [df for df in series_list if isinstance(df, pl.LazyFrame)]
    if lazy_dfs:
        df = pl.concat(lazy_dfs)
        df_normalized = normalize_multimapping_umitools_group_output(df)
        count_group = pivot_normalized_count(df_normalized, dedup_umi=dedup_umi)
    else:
        count_group = None
    matrix = (
        pd.concat(
            [i for i in series_list if isinstance(i, pd.Series)]
            + ([count_group] if count_group is not None else []),
            axis=0,
        )
        .fillna(0)
        .astype(int)
    ).loc[samples]
    return matrix


def normalize_multimapping_umitools_group_output(df: pl.LazyFrame) -> pl.LazyFrame:
    # 1. Build per-read dataframe and construct gene co-occurrence graph
    df_per_read = consolidate_read_id(df)
    # all_genes = (
    #     df_per_read.select(pl.col("gene").explode())
    #     .unique()
    #     .collect()["gene"]
    #     .to_list()
    # )
    # Build graph of overlapping features
    edges_df = per_read_to_edges(df_per_read)
    graph = edges_to_graph(edges_df)

    # 2. Create component names using simplify_gene_list
    gene_to_component = {
        gene: "|".join(sorted(component))
        for component in nx.connected_components(graph)
        for gene in component
    }
    # Check that each read_id has one normalized_gene
    assert (
        df_per_read.select(
            pl.col("gene").implode().replace(gene_to_component).n_unique()
        )
        .max()
        .collect()
        .item()
        == 1
    ), "Each read_id should have exactly one normalized gene name"

    df_normalized = df_per_read.with_columns(
        pl.col("gene").list.first().replace(gene_to_component)
    ).unique(["sample", "read_id"])
    return df_normalized


def pivot_normalized_count(df: pl.LazyFrame, dedup_umi: bool = True) -> pd.DataFrame:
    return (
        df.group_by(["sample", "gene"])
        .agg(pl.n_unique("final_umi") if dedup_umi else pl.len())
        .select(["sample", "gene", "final_umi"])
        .collect()
        .pivot(index="sample", on="gene", values="final_umi")
        .fill_null(0)
        .to_pandas()
        .set_index("sample")
    )


def id2symbol_from_gencode_gtf(gtf: str) -> dict[str, str]:
    id2symbol: dict[str, str] = {}
    gtf_df = pl.read_csv(
        gtf,
        comment_prefix="#",
        has_header=False,
        new_columns=[
            "seqname",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "frame",
            "attribute",
        ],
        separator="\t",
    ).to_pandas()
    gtf_df["gene_id"] = gtf_df.attribute.str.extract(r'gene_id "(.*?)";')
    gtf_df["gene_symbol"] = gtf_df.attribute.str.extract(r'gene_name "(.*?)";')
    for gene_id, gene_symbol in gtf_df[["gene_id", "gene_symbol"]].values:
        if not gene_id:
            continue
        gene_id_nover = gene_id.split(".")[0]
        if gene_id_nover in id2symbol:
            assert id2symbol[gene_id_nover] == gene_symbol
        else:
            id2symbol[gene_id_nover] = gene_symbol
    return id2symbol


def read_counts(
    input_: str | list[str], gtf: str | None = None, dedup_umi: bool = True
) -> pd.DataFrame:
    """
    Read one or more count files into a sample x gene count matrix.

    Args:
        input_: A single file path, a glob pattern, or a list of file paths.
        gtf: Optional path to a GTF file for gene_id→gene_symbol mapping.
             If provided, gene IDs (without version) will be renamed to symbols.

    Returns:
        pd.DataFrame: Rows are samples, columns are genes (symbols if gtf provided),
                      values are integer UMI counts.
    """
    files = expand_path_to_list(input_)
    # process each file
    matrix = process_count_file(files, dedup_umi=dedup_umi)
    # combine into matrix: samples × genes

    # rename columns to gene symbols if map provided
    if gtf:
        id2symbol = id2symbol_from_gencode_gtf(gtf)
        matrix.columns = [
            "|".join(
                simplify_gene_list(
                    [id2symbol.get(i.split(".")[0], i) for i in col.split("|")]
                )
            )
            for col in matrix.columns
        ]
        # Aggregate by gene symbol if multiple IDs map to the same symbol
        matrix = (
            matrix.transpose()
            .groupby(level=0, sort=False)
            .mean()
            .transpose()
            .astype(int)
        )
    # reorder columns by average rel ab
    rel_ab = matrix.divide(matrix.sum(axis=1), axis=0)
    matrix = matrix[rel_ab.mean(axis=0).sort_values(ascending=False).index]

    return matrix
