#!/usr/bin/env python3
"""
Extract region sequences for the final selected primer set ("product" stage).

Reads:
  - Selected primer pairs TSV (output from the SADDLE stage).
  - Primer metadata TSV (output from the subset stage).
  - Reference transcript FASTA (packaged, resolved by --species; or user-specified).

Writes:
  - FASTA of amplified regions for each selected primer pair.

Typical columns:
  - Selected TSV: transcriptID, primerIndex, sequenceLeft, sequenceRight
  - Primer info TSV (metadata): includes amplicon_index, geneID, geneSymbol,
    transcriptID, L_seq, L_start, R_seq, R_start, amplicon_size
"""

import gzip
from pathlib import Path
from typing import Literal

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .utils import resolve_reference_path


def _open_fasta_auto(path: Path):
    """Return a readable text handle for FASTA (gz or plain).

    Args:
        path: Path to the FASTA file (may be .gz).

    Returns:
        A file-like text handle.
    """
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def load_reference_transcripts(ref_file: Path) -> dict[str, str]:
    """Load reference transcript sequences from a FASTA.

    Args:
        ref_file: Path to transcript FASTA (plain or gzipped).

    Returns:
        Mapping transcript_id -> sequence string.
    """
    with _open_fasta_auto(ref_file) as handle:
        return {rec.id: str(rec.seq) for rec in SeqIO.parse(handle, "fasta")}


def load_primer_info(primer_file: Path) -> pd.DataFrame:
    """Load and normalize primer metadata.

    Args:
        primer_file: TSV with primer metadata from the subset stage.

    Returns:
        DataFrame indexed by amplicon_index with normalized columns:
        [gene_id, gene_symbol, transcript_id, l_seq, l_start, r_seq, r_start, amplicon_size]
    """
    df = pd.read_csv(primer_file, sep="\t")
    cols_needed = [
        "amplicon_index",
        "geneID",
        "geneSymbol",
        "transcriptID",
        "L_seq",
        "L_start",
        "R_seq",
        "R_start",
        "amplicon_size",
    ]
    missing = [c for c in cols_needed if c not in df.columns]
    if missing:
        raise ValueError(f"Primer info missing required columns: {missing}")

    primer_df = df[cols_needed].copy()
    primer_df.columns = [
        "amplicon_index",
        "gene_id",
        "gene_symbol",
        "transcript_id",
        "l_seq",
        "l_start",
        "r_seq",
        "r_start",
        "amplicon_size",
    ]
    return primer_df.set_index("amplicon_index")


def _coerce_selected_columns(selected_df: pd.DataFrame) -> pd.DataFrame:
    """Coerce/rename columns from selected TSV.

    Args:
        selected_df: DataFrame read from the SADDLE output TSV.

    Returns:
        DataFrame with columns: transcript_id, primer_index, l_seq, r_seq
    """
    # Accept either exact names or position-based fallback.
    cols_needed = ["transcriptID", "primerIndex", "sequenceLeft", "sequenceRight"]
    cols_renamed = ["transcript_id", "primer_index", "l_seq", "r_seq"]
    if set(cols_needed).issubset(selected_df.columns):
        out = selected_df[cols_needed].copy()
    else:
        # Fallback to positional rename (4 columns expected).
        out = selected_df.copy()
        if len(out.columns) < 4:
            raise ValueError("Selected TSV must have at least 4 columns.")
        out = out.iloc[:, :4]
        out.columns = cols_needed

    out.columns = cols_renamed
    return out


def product(
    selected_tsv: Path,
    primer_info_tsv: Path,
    output_fasta: Path,
    species: Literal["mouse", "human"] | None = None,
    reference: Path | None = None,
) -> int:
    """Extract and save selected regions to a FASTA file.

    Args:
        selected_tsv: TSV from the SADDLE stage (selected primer pairs). Output by foli-primer saddle.
        primer_info_tsv: Primer metadata TSV from the subset stage. Output by foli-primer subset.
        output_fasta: Output FASTA path to write regions.
        reference_fasta: Transcript FASTA (plain or gzipped).

    Returns:
        Exit code 0 on success.
    """
    selected_df = _coerce_selected_columns(pd.read_csv(selected_tsv, sep="\t"))

    primer_df = load_primer_info(primer_info_tsv)

    if species is not None:
        ref_path = resolve_reference_path(None, species)
    elif reference is not None:
        ref_path = resolve_reference_path(reference, None)
    else:
        raise ValueError("Either species or reference must be provided.")
    transcript_pool = load_reference_transcripts(ref_path)

    success = 0
    errors = 0
    records: list[SeqRecord] = []

    for _, row in selected_df.iterrows():
        tx_id = str(row["transcript_id"])
        primer_idx = row["primer_index"]
        l_seq = row["l_seq"]
        r_seq = row["r_seq"]

        try:
            info = primer_df.loc[primer_idx]
        except KeyError:
            print(
                f"Error [{tx_id}]: missing primer_index in primer info -> {primer_idx}"
            )
            errors += 1
            continue

        cand_l_seq = info["l_seq"]
        cand_r_seq = info["r_seq"]
        cand_l_start = int(info["l_start"])  # 1-based
        cand_r_start = int(info["r_start"])  # 1-based
        amplicon_size = int(info["amplicon_size"])

        if not (cand_l_seq == l_seq and cand_r_seq == r_seq):
            print(f"Error [{tx_id}]: primer sequence mismatch")
            errors += 1
            continue

        try:
            seq_full = transcript_pool[tx_id]
        except KeyError:
            print(f"Error [{tx_id}]: transcript not found in reference FASTA")
            errors += 1
            continue

        # Convert to 0-based inclusive slice [L_start-1 : R_start)
        region = seq_full[(cand_l_start - 1) : cand_r_start]
        if len(region) != amplicon_size:
            print(
                f"Error [{tx_id}]: length mismatch (expected {amplicon_size}, got {len(region)})"
            )
            errors += 1
            continue

        rec = SeqRecord(
            Seq(region),
            id=str(primer_idx),
            description=f"gene:{info['gene_symbol']} transcript:{tx_id}",
        )
        records.append(rec)
        success += 1
        print(f"Success [{tx_id}]: {len(region)} bp")

    if records:
        output_fasta.parent.mkdir(parents=True, exist_ok=True)
        SeqIO.write(records, output_fasta, "fasta")

    print(f"\nSummary: {success} successful extractions, {errors} errors")
    print(f"Output written to: {output_fasta}")
    return 0
