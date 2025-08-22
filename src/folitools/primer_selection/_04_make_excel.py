import pandas as pd

from .utils import get_prefixes


def _generate_idt_order_file(
    res_df: pd.DataFrame,
    output_idt_order: str,
    idt_pool_prefix: str,
) -> None:
    """
    Generate IDT ordering Excel file with pool assignments.

    Args:
        res_df: DataFrame with primer information sorted by Group, geneSymbol, etc.
                Must contain columns: primer_sequence_to_order_forward, primer_sequence_to_order_reverse
        output_idt_order: Path to output IDT ordering Excel file.
        idt_pool_prefix: Prefix for IDT pool names.
    """
    # Create list of primers in order: forward first, then reverse for each gene
    idt_primers = []

    for _, row in res_df.iterrows():
        # Add forward primer
        idt_primers.append(
            {
                "Sequence": row["primer_sequence_to_order_forward"],
                "Gene": row["geneSymbol"],
                "Group": row["Group"],
                "Direction": "forward",
            }
        )

        # Add reverse primer
        idt_primers.append(
            {
                "Sequence": row["primer_sequence_to_order_reverse"],
                "Gene": row["geneSymbol"],
                "Group": row["Group"],
                "Direction": "reverse",
            }
        )

    # Create DataFrame
    idt_df = pd.DataFrame(idt_primers)

    # Assign pool names based on 384 primer limit
    max_primers_per_pool = 384
    total_primers = len(idt_df)

    pool_names = []
    for i in range(total_primers):
        pool_number = i // max_primers_per_pool

        if total_primers <= max_primers_per_pool:
            # Single pool, no suffix
            pool_name = idt_pool_prefix
        else:
            # Multiple pools, add suffix
            pool_name = f"{idt_pool_prefix}-p{pool_number + 1}"

        pool_names.append(pool_name)

    idt_df["Pool name"] = pool_names

    # Create final output with just the two required columns
    idt_output = idt_df[["Pool name", "Sequence"]]

    # Save to Excel
    idt_output.to_excel(output_idt_order, index=False)

    # Print summary
    pool_counts = idt_df["Pool name"].value_counts().sort_index()
    print(f"IDT ordering file created with {len(pool_counts)} pool(s):")
    for pool, count in pool_counts.items():
        print(f"  {pool}: {count} primers")


def summary(
    input_: str,
    primer_selection: str,
    primer_info: str,
    output: str,
    has_linker: bool = False,
    output_idt_order: str | None = None,
    idt_pool_prefix: str = "pool",
) -> pd.DataFrame:
    """
    Generate an Excel file for primer ordering.
    Args:
        input_ (str): Path to gene info tsv file.
        primer_selection (str): Path to selected primers tsv file.
        primer_info (str): Path to candidate primer info tsv file.
        output (str): Path to output Excel file (should end with .xlsx).
        has_linker (bool): Whether to include linker sequences in primers.
        output_idt_order (str | None): Optional path to output IDT ordering Excel file.
        idt_pool_prefix (str): Prefix for IDT pool names (default: "pool").
    """
    selection = pd.read_table(primer_selection)
    primer_info_df = pd.read_table(primer_info)
    df_gene = pd.read_table(input_).rename(
        {"group": "Group", "gene": "geneSymbol"}, axis=1
    )

    # primer_info["amplicon_index"].is_unique, selection["primerIndex"].is_unique

    assert df_gene["geneSymbol"].is_unique, "geneSymbol should be unique"

    res_df = pd.merge(
        selection,
        primer_info_df,
        left_on="primerIndex",
        right_on="amplicon_index",
        how="left",
    ).merge(df_gene, on="geneSymbol", how="left")
    res_df["Chosen Index"] = res_df.primerIndex.map(lambda x: int(x.split("_")[1]))
    # reorder rows
    # res_df = (
    #     res_df.set_index("geneSymbol")
    #     .loc[[i for i in gene_names if i in res_df.geneSymbol.tolist()]]
    #     .reset_index()
    # )
    # take columns
    res_df = res_df[
        [
            "Group",
            "geneSymbol",
            "geneID",
            "Chosen Index",
            "amplicon_index",
            "L_seq",
            "R_seq",
            "pool",
        ]
    ]

    # format to excel like:
    # amplicon_index  L_seq  R_seq  primer_sequence_to_order_forward  primer_sequence_to_order_reverse
    # ENSMUST00000001184_0  TTTTCACCTACCCCTGCTGTGTTCG  GCAGGGAGACAAAACACTGTAGCCA  CCTACACGACGCTCTTCCGATCTNNNNNNTTTTCACCTACCCCTGCTGTGTTCG  TCAGACGTGTGCTCTTCCGATCTNNNNNNGCAGGGAGACAAAACACTGTAGCCA
    # ENSMUST00000001184_1  CATCAGGTCTGTTCACTGCACCCAA  ACAGTGACTCAGTCCAGCTTCCAGA  CCTACACGACGCTCTTCCGATCTNNNNNNCATCAGGTCTGTTCACTGCACCCAA  TCAGACGTGTGCTCTTCCGATCTNNNNNNACAGTGACTCAGTCCAGCTTCCAGA
    FWD_PREFIX, REV_PREFIX = get_prefixes(has_linker)

    res_df["primer_sequence_to_order_forward"] = FWD_PREFIX + res_df["L_seq"]
    res_df["primer_sequence_to_order_reverse"] = REV_PREFIX + res_df["R_seq"]
    res_df = res_df.sort_values(by=["pool", "Group", "geneSymbol", "Chosen Index"])
    res_df.to_excel(output, index=False)
    print(res_df["Group"].value_counts())

    # Generate IDT ordering file if requested
    if output_idt_order is not None:
        _generate_idt_order_file(res_df, output_idt_order, idt_pool_prefix)

    return res_df
