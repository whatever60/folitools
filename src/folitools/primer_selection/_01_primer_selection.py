#!/usr/bin/env python3
"""
Primer Selection (subset logic)

This module implements the filtering logic to subset primer sequences and
their metadata given a species, an amplicon size range, and a user-supplied
gene table.

Input data discovery
--------------------
Packaged inputs are resolved from:
    folitools/primer_selection/data/output_<species>/

Expected Parquet files (UTF-8; Parquet; schema noted below):
    1) 01.primer3_pass.<min>_<max>.primer_info.parquet
       Columns (minimum):
         - transcript_index : int or str (join key to primer seq)
         - geneID           : str
         - geneSymbol       : str
         - amplicon_index   : int or str (unique per amplicon)
         - L_seq            : str (forward primer)
         - R_seq            : str (reverse primer)
         - chrN             : str (kept as string dtype if present)
       Additional columns are preserved.

    2) 02.SADDLE_input.<min>_<max>.primer_0.parquet
       Columns (minimum):
         - <col0> transcriptID     : str (matches transcript metadata)
         - <col1> transcript_index : int or str (join key to primer info)
         - primerIndex             : int or str (maps to amplicon_index)
       Additional columns are preserved.

    3) 02.SADDLE_input.<min>_<max>.transcript_metadata.parquet
       Columns (minimum; also accessed by position):
         - geneID        : str
         - geneSymbol    : str  (assumed at positional index 2)
         - transcriptID  : str  (assumed at positional index 3)
       The module uses positional access for geneSymbol (col 2) and
       transcriptID (col 3); ensure these align with your files.

User-supplied TSV
-----------------
Gene table TSV (required by the CLI):
  Required columns:
    - gene   : gene symbol (unique preferred)
    - group  : arbitrary grouping/annotation
  Optional columns:
    - primer_fwd : explicit forward primer to enforce
    - primer_rev : explicit reverse primer to enforce

Outputs (TSV)
-------------
Written to --output-dir (defaults to "."):
  - candidate_primer.<min>_<max>.primer_sequence.tsv
  - candidate_primer.<min>_<max>.primer_info.tsv
  - candidate_primer.<min>_<max>.transcript_metadata.tsv

Note:
This module exposes a callable `subset(...)`. The CLI entrypoint
lives in `folitools/primer_selection/cli.py` (subcommand: `subset`).
"""

from pathlib import Path

import pandas as pd
import polars as pl

from .data_dir_resolve import _data_dir_for_species


def read_gene_metadata(file_path: Path) -> pd.DataFrame:
    """Read the gene table TSV with optional explicit primer pairs.

    Args:
        file_path: Path to a TSV containing at least columns ``gene`` and ``group``.
            Optional columns: ``primer_fwd``, ``primer_rev``.

    Returns:
        A pandas DataFrame of the gene table.

    Raises:
        ValueError: If required columns are missing.
    """
    df = pd.read_csv(file_path, sep="\t")
    required = {"gene", "group"}
    if not required.issubset(df.columns):
        raise ValueError(f"Gene table must contain columns: {required}")
    return df


def _load_and_filter_data(
    transcript_df: pd.DataFrame,
    primer_seq_df: pd.DataFrame,
    primer_info_df: pd.DataFrame,
    target_genes: set[str],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Stage-wise filtering using gene symbols and join keys.

    Filtering steps:
      1) Keep transcripts whose gene symbol (positional column index 2) is in target_genes.
      2) Keep primer-seq rows whose transcriptID (primer-seq col 0) is among the
         retained transcripts' transcriptID (positional index 3).
      3) Keep primer-info rows whose transcript_index (primer-info col 0)
         is among the primer-seq transcript_index (primer-seq col 1).
      4) Constrain transcripts to those whose geneID intersects primer-info geneID.
         Warn if selected transcript identifiers are non-unique.

    Args:
        transcript_df: Transcript metadata table.
        primer_seq_df: Primer sequence table.
        primer_info_df: Primer metadata table.
        target_genes: Target gene symbols.

    Returns:
        (filtered_transcripts, filtered_primer_seq, filtered_primer_info)
    """
    print(f"Loaded {len(transcript_df)} transcript records")
    gene_col = transcript_df.columns[2]
    txid_col = transcript_df.columns[3]

    transcripts_keep = transcript_df[transcript_df[gene_col].isin(target_genes)]
    print(
        f"Found {len(transcripts_keep)} transcripts for {len(target_genes)} target genes"
    )

    # Primer sequences by transcriptID (col 0) and carry transcript_index (col 1)
    seq_txid_col = primer_seq_df.columns[0]
    seq_txidx_col = primer_seq_df.columns[1]

    primer_seq_keep = primer_seq_df[
        primer_seq_df[seq_txid_col].isin(set(transcripts_keep[txid_col]))
    ]
    print(f"Found {len(primer_seq_keep)} primer sequences for target transcripts")

    # Primer info by transcript_index (col 0)
    info_txidx_col = primer_info_df.columns[0]
    primer_info_keep = primer_info_df[
        primer_info_df[info_txidx_col].isin(set(primer_seq_keep[seq_txidx_col]))
    ]
    print(f"Found {len(primer_info_keep)} primer metadata records for target indices")

    # Sanity: intersect transcripts by geneID with primer-info geneID
    if "geneID" in transcripts_keep.columns and "geneID" in primer_info_keep.columns:
        transcripts_keep = transcripts_keep[
            transcripts_keep["geneID"].isin(primer_info_keep["geneID"])
        ]

    # Uniqueness warnings
    if not (
        transcripts_keep["geneID"].is_unique
        and transcripts_keep["geneSymbol"].is_unique
        and transcripts_keep["transcriptID"].is_unique
    ):
        print(
            "Warning: filtered transcript identifiers are not unique within the selection."
        )

    return transcripts_keep, primer_seq_keep, primer_info_keep


def _post_filter_by_explicit_primers(
    primer_info_df: pd.DataFrame,
    primer_seq_df: pd.DataFrame,
    gene_meta_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Optionally restrict to explicit primer pairs per gene.

    Args:
        primer_info_df: Filtered primer metadata, requires columns:
            ``geneSymbol``, ``L_seq``, ``R_seq``, ``amplicon_index``.
        primer_seq_df: Filtered primer sequences, requires column ``primerIndex``.
        gene_meta_df: Gene table; rows may contain ``primer_fwd`` and/or ``primer_rev``.

    Returns:
        (primer_info_df_filtered, primer_seq_df_filtered)

    Raises:
        ValueError: If an explicit (gene, fwd, rev) combination is not found.
    """
    keep_idx: set = set()

    for _, row in gene_meta_df.iterrows():
        fwd = row.get("primer_fwd")
        rev = row.get("primer_rev")
        matches = primer_info_df.query("geneSymbol == @row['gene']")

        if pd.notna(fwd):
            matches = matches[matches["L_seq"] == fwd]
        if pd.notna(rev):
            matches = matches[matches["R_seq"] == rev]
        if matches.empty:
            raise ValueError(
                f"Primer pair not found for gene {row['gene']} with fwd={fwd}, rev={rev}"
            )
        keep_idx.update(matches["amplicon_index"])

    info_f = primer_info_df[primer_info_df["amplicon_index"].isin(keep_idx)]
    seq_f = primer_seq_df[primer_seq_df["primerIndex"].isin(keep_idx)]
    return info_f, seq_f


def _validate_output_file_extensions(*file_paths: Path) -> None:
    """Validate that output file paths have acceptable extensions.

    Args:
        *file_paths: Variable number of Path objects to validate.

    Raises:
        ValueError: If any file path doesn't end with an acceptable extension.
    """
    valid_extensions = {".tsv", ".txt", ".tsv.gz", ".txt.gz"}

    for file_path in file_paths:
        # Handle compressed files (.tsv.gz, .txt.gz)
        if len(file_path.suffixes) >= 2 and file_path.suffixes[-2] in {".tsv", ".txt"}:
            ext = "".join(file_path.suffixes[-2:])
        else:
            ext = file_path.suffix

        if ext not in valid_extensions:
            raise ValueError(
                f"Output file '{file_path}' must end with one of {valid_extensions}, got: {ext}"
            )


def subset(
    gene_table_file: Path,
    species: str,
    amplicon_size_range: tuple[int, int],
    output_primer_sequence: Path,
    output_primer_info: Path,
) -> int:
    """Run the subset operation end-to-end.

    Args:
        species: Species name driving packaged input discovery (e.g., "mouse").
        amplicon_size_range: Range like "320-380" (also accepts "320:380").
        gene_table_file: Path to the user TSV with columns ``gene`` and ``group``,
            optionally ``primer_fwd`` and ``primer_rev``.
        output_primer_sequence: Path to write primer sequence TSV (must end with .tsv/.txt/.tsv.gz/.txt.gz).
        output_primer_info: Path to write primer info TSV (must end with .tsv/.txt/.tsv.gz/.txt.gz).

    Returns:
        Exit code: 0 on success.

    Raises:
        FileNotFoundError: If required packaged Parquet files are missing.
        ValueError: On invalid inputs or explicit primer mismatches.
    """
    # amin, amax = _parse_amplicon_size_range(amplicon_size_range)
    amin, amax = amplicon_size_range

    # Validate output file extensions
    _validate_output_file_extensions(output_primer_sequence, output_primer_info)

    # Resolve packaged inputs
    suffix = f"{amin}_{amax}"
    data_dir = _data_dir_for_species(species)
    p_info = data_dir / f"01.primer3_pass.{suffix}.primer_info.parquet"
    p_seq = data_dir / f"02.SADDLE_input.{suffix}.primer_0.parquet"
    p_tx = data_dir / f"02.SADDLE_input.{suffix}.transcript_metadata.parquet"
    for p in (p_info, p_seq, p_tx):
        if not p.exists():
            raise FileNotFoundError(f"Packaged input not found: {p}")

    # Ensure output directories exist
    output_primer_sequence.parent.mkdir(parents=True, exist_ok=True)
    output_primer_info.parent.mkdir(parents=True, exist_ok=True)

    # Load inputs
    transcript_df = pl.read_parquet(p_tx).to_pandas()
    primer_seq_df = pl.read_parquet(p_seq).to_pandas()
    primer_info_df = pl.read_parquet(p_info).to_pandas()

    # Ensure chrN is string-like, if present
    if "chrN" in primer_info_df.columns:
        primer_info_df["chrN"] = primer_info_df["chrN"].astype("string")

    # Load gene table
    gene_meta_df = read_gene_metadata(gene_table_file)
    if not gene_meta_df["gene"].is_unique:
        print("Warning: Gene table contains non-unique gene symbols.")
    target_genes = set(gene_meta_df["gene"])
    print(f"Loaded {len(target_genes)} target gene symbols")

    # Stage-wise filtering
    tx_f, seq_f, info_f = _load_and_filter_data(
        transcript_df=transcript_df,
        primer_seq_df=primer_seq_df,
        primer_info_df=primer_info_df,
        target_genes=target_genes,
    )

    # Optional explicit pairs
    genes_with_any = gene_meta_df[
        gene_meta_df["gene"].isin(info_f["geneSymbol"].unique())
    ]
    if {"primer_fwd", "primer_rev"}.intersection(genes_with_any.columns):
        # Only restrict if at least one explicit value is present in *any* row
        any_explicit = (
            genes_with_any[["primer_fwd", "primer_rev"]].notna().any(axis=None)
        )
        if any_explicit:
            info_f, seq_f = _post_filter_by_explicit_primers(
                info_f, seq_f, genes_with_any
            )

    # Save outputs
    seq_f.to_csv(output_primer_sequence, sep="\t", index=False)
    info_f.to_csv(output_primer_info, sep="\t", index=False)
    # tx_f.to_csv(out_tx, sep="\t", index=False)

    print(f"Wrote: {output_primer_sequence}")
    print(f"Wrote: {output_primer_info}")
    # print(f"Wrote: {out_tx}")
    print("Subset complete.")
    return 0
