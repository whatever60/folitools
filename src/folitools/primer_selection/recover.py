"""Recover primer FASTA files from IDT order Excel."""

import logging
import os
import shutil
import subprocess
from io import StringIO
from pathlib import Path
from typing import Literal

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Local imports
from .recover_plot import make_report
from .utils import get_prefixes, resolve_reference_path

__all__ = ["recover"]

logger = logging.getLogger(__name__)


def _read_idt_excel(order_excel: Path, robust: bool = True) -> pd.DataFrame:
    """Read the IDT order Excel file and extract two columns: pool and sequence.

    The function is tolerant to header variations and discards obvious non-sequence rows.

    Args:
        order_excel: Path to the IDT Excel file.

    Returns:
        DataFrame with columns ["pool", "sequence"].
    """
    if not robust:
        return pd.read_excel(order_excel, usecols=[0, 1], names=["pool", "sequence"])
    df_raw = pd.read_excel(order_excel, usecols=[0, 1], header=0)
    df_raw = df_raw.rename(
        columns={df_raw.columns[0]: "pool", df_raw.columns[1]: "sequence"}
    )
    # Drop non-string and apparent headers/NaNs
    df_raw = df_raw.dropna(subset=["sequence"])
    df_raw = df_raw[df_raw["sequence"].map(lambda x: isinstance(x, str))]
    df_raw["sequence"] = df_raw["sequence"].str.strip().str.upper()
    # Keep rows that look like DNA with N's
    mask_dna = df_raw["sequence"].str.fullmatch(r"[ACGTN]+", na=False)
    df = df_raw[mask_dna].copy()
    if df.empty:
        raise ValueError(
            "No DNA-like sequences found in the first two columns of the Excel file."
        )
    return df.reset_index(drop=True)


def _classify_primers(
    df_idt: pd.DataFrame, fwd_prefix: str, rev_prefix: str
) -> pd.DataFrame:
    """Classify each row as forward or reverse primer and extract gene-specific sequences.

    Args:
        df_idt: DataFrame with columns ["pool", "sequence"].
        fwd_prefix: Forward primer universal prefix (includes NNNNNN and optional linker).
        rev_prefix: Reverse primer universal prefix (includes NNNNNN and optional linker).

    Returns:
        DataFrame with columns: [pool, primer_type, primer_seq_full, primer_seq].

    Raises:
        ValueError: If any sequence doesn't match either prefix.
    """
    primers = []
    unmatched: list[str] = []

    for _, row in df_idt.iterrows():
        seq = row["sequence"].upper()
        if seq.startswith(fwd_prefix):
            gene_part = seq[len(fwd_prefix) :]
            primers.append(
                {
                    "pool": row["pool"],
                    "primer_type": "fwd",
                    "primer_seq_full": seq,
                    "primer_seq": gene_part,
                }
            )
        elif seq.startswith(rev_prefix):
            gene_part = seq[len(rev_prefix) :]
            primers.append(
                {
                    "pool": row["pool"],
                    "primer_type": "rev",
                    "primer_seq_full": seq,
                    "primer_seq": gene_part,
                }
            )
        else:
            unmatched.append(str(row.get("pool", "<NA>")))

    if unmatched:
        raise ValueError(
            "Some primers do not match the expected prefixes. Offending pools: "
            + ", ".join(unmatched)
        )

    primer_df = pd.DataFrame(primers)
    # Sanity: primer gene-specific sequences should be unique
    dup_counts = primer_df["primer_seq"].value_counts()
    dups = dup_counts[dup_counts > 1]
    if not dups.empty:
        logger.warning(
            "Duplicate primer gene sequences detected: %s",
            ", ".join(dups.index.tolist()),
        )
    return primer_df


def _ensure_seqkit_available() -> None:
    """Raise if `seqkit` is not on PATH."""
    if shutil.which("seqkit") is None:
        raise RuntimeError(
            "`seqkit` not found on PATH. Please install seqkit and try again."
        )


def _run_seqkit_locate(
    primer_seqs: list[str], ref_fasta: Path, threads: int = 1
) -> pd.DataFrame:
    """Run `seqkit locate` for a batch of patterns and return its table as a DataFrame.

    Args:
        primer_seqs: List of gene-specific primer sequences (A/C/G/T/N, case-insensitive).
        ref_fasta: Path to the transcriptome FASTA file.
        threads: Number of threads for seqkit; defaults to `os.cpu_count()`.

    Returns:
        A DataFrame parsed from `seqkit locate` TSV output.

    Raises:
        RuntimeError: If `seqkit locate` fails.
    """
    _ensure_seqkit_available()
    t = threads or (os.cpu_count() or 1)

    args: list[str] = [
        "seqkit",
        "locate",
        "--ignore-case",
        "--max-mismatch",
        "0",
        "--use-fmi",
    ]
    for p in primer_seqs:
        args.extend(["--pattern", p])
    args.extend(["--threads", str(t), str(ref_fasta)])

    logger.info("Running: %s", " ".join(args))
    proc = subprocess.run(args, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"seqkit locate failed: {proc.stderr}")

    df = pd.read_csv(StringIO(proc.stdout), sep="\t")

    # Normalize seqID column name across seqkit versions
    if "seqID" not in df.columns:
        raise RuntimeError(
            "Unexpected seqkit output: missing 'seqID' or 'sequence' column."
        )

    return df


def _enrich_locate_df(
    locate_df: pd.DataFrame, primer_info: pd.DataFrame
) -> pd.DataFrame:
    """Join locate results with primer metadata; parse transcript/gene fields from FASTA headers.

    Expects FASTA headers in the format compatible with the user's packaged data, where
    "seqID" is a pipe-separated string and indices [0] = transcript_id, [1] = gene_id,
    [5] = gene_symbol.
    """
    df = locate_df.copy()
    df["primer_seq"] = df["pattern"].str.upper()

    # Parse pipe-separated identifiers safely
    parts = df["seqID"].astype(str).str.split("|", expand=True)
    df["transcript_id"] = parts[0]
    df["gene_id"] = parts[1].fillna("")
    df["gene_symbol"] = parts[5] if parts.shape[1] > 5 else ""

    merged = pd.merge(df, primer_info, how="left", on="primer_seq")
    missing = merged[merged["pool"].isna()]
    if not missing.empty:
        raise RuntimeError(
            "Internal error: locate patterns didn't map back to primers."
        )

    locate_df_final = merged[
        [
            "primer_seq",
            "primer_seq_full",
            "primer_type",
            "transcript_id",
            "gene_id",
            "gene_symbol",
            "start",
            "end",
            "strand",
            "pool",
        ]
    ].copy()

    # Sanity check: Ideally, no duplicate primer_seq and transcript_id pairs. In reality,
    # a primer sequence can bind to multiple locations on the same transcript if the
    # transcript sequence is repetitive.
    for (primer_seq, transcript_id), group in locate_df_final.groupby(
        ["primer_seq", "transcript_id"]
    ):
        if len(group) > 1:
            logger.info(
                f"Duplicate match found for primer sequence {primer_seq} and transcript ID "
                f"{transcript_id}: {len(group)} matches."
            )
    # Sanity check: Ideally, a primer_seq matches to one gene_id only. In reality, a
    # primer sequence could bind to multiple locations on multiple transcripts of the same gene.
    for primer_seq, group in locate_df_final.groupby("primer_seq"):
        gene_ids = group["gene_id"].unique()
        if len(gene_ids) > 1:
            logger.info(
                f"Primer sequence {primer_seq} matches to {len(gene_ids)} different gene "
                f"IDs: {', '.join(gene_ids)}."
            )
    return merged


def _build_amplicons(locate_final: pd.DataFrame) -> pd.DataFrame:
    """Enumerate amplicons by pairing fwd/rev hits on each transcript.

    The pairing logic is as follows:
      - Forward primer should be equivalent to transcript sequence; reverse primer RC to transcript.
      - Opposite strands required.
      - Upstream end must be before downstream start on the transcript coordinate.

    Args:
        locate_final: Joined table with primer/transcript hit info.

    Returns:
        DataFrame with amplicon information.
    """
    rows: list[dict] = []

    for tid, group in locate_final.groupby("transcript_id"):
        g_fwd = group[group["primer_type"] == "fwd"]
        g_rev = group[group["primer_type"] == "rev"]
        # if g_fwd.empty or g_rev.empty:
        #     continue
        for _, rf in g_fwd.iterrows():
            for _, rr in g_rev.iterrows():
                flipped = False
                up, down = rf, rr

                # Require opposite strands
                if rf["strand"] == rr["strand"]:
                    continue

                # Cases by strand orientation as in the notebook
                if rf["strand"] == "-" and rr["strand"] == "+":
                    if rf["start"] > rr["end"]:
                        up, down = rr, rf
                        flipped = True
                    else:
                        continue
                else:
                    if not (rf["end"] < rr["start"]):
                        continue
                # NOTE: After the above check, we made it such that:
                # - Forward primers are always equivalent to the transcript sequences.
                # - Reverse primers are always reverse-complemented to the transcript sequences.
                amplicon = {
                    "primer_seq_fwd": rf["primer_seq"],
                    "primer_seq_rev": rr["primer_seq"],
                    "transcript_id": tid,
                    "gene_id": up["gene_id"],
                    "gene_symbol": up["gene_symbol"],
                    "start_up": int(up["start"]),  # 1-based, inclusive
                    "end_up": int(up["end"]),
                    "start_down": int(down["start"]),
                    "end_down": int(down["end"]),
                    "pool_fwd": rf["pool"],
                    "pool_rev": rr["pool"],
                    "flipped": flipped,
                }
                amplicon["amplicon_length"] = (
                    amplicon["end_down"] - amplicon["start_up"] + 1
                )
                rows.append(amplicon)

    amplicon_all = pd.DataFrame(rows)
    if amplicon_all.empty:
        raise ValueError("No amplicons were formed. Check inputs and prefixes.")
    return amplicon_all


def _group_pairs(amplicon_sub: pd.DataFrame) -> pd.DataFrame:
    """Group by primer pair and summarize mappings across transcripts/genes.

    Returns a DataFrame with one row per (forward, reverse) pair with stable fields
    and a semicolon-joined "multi_mapping" column.
    """
    if amplicon_sub.empty:
        return pd.DataFrame(
            columns=[
                "primer_seq_fwd",
                "primer_seq_rev",
                "transcript_id",
                "gene_id",
                "gene_symbol",
                "start_up",
                "end_up",
                "start_down",
                "end_down",
                "pool_fwd",
                "pool_rev",
                "num_transcripts",
                "num_genes",
                "transcript_id_all",
            ]
        )

    rows: list[dict] = []
    grouped = amplicon_sub.groupby(["primer_seq_fwd", "primer_seq_rev"], sort=False)
    for (fwd, rev), g in grouped:
        pool_fwd_vals = g["pool_fwd"].unique()
        pool_rev_vals = g["pool_rev"].unique()
        assert len(pool_fwd_vals) == 1, f"Multiple pool_fwd for primers: {fwd}, {rev}"
        assert len(pool_rev_vals) == 1, f"Multiple pool_rev for primers: {fwd}, {rev}"
        # if len(pool_fwd_vals) != 1 or len(pool_rev_vals) != 1:
        #     logger.warning(
        #         "Multiple pools found for primer pair %s / %s: L=%s, R=%s",
        #         fwd,
        #         rev,
        #         pool_fwd_vals.tolist(),
        #         pool_rev_vals.tolist(),
        #     )

        g_sorted = g.sort_values(by=["gene_id", "transcript_id"])
        parts = []
        for _, r in g_sorted.iterrows():
            part = (
                f"{r.transcript_id}|{r.gene_id}|{r.gene_symbol}|"
                f"{r.start_up}|{r.end_up}|{r.start_down}|{r.end_down}|"
                f"{r.amplicon_length}"
            )
            parts.append(part)

        rows.append(
            {
                "primer_seq_fwd": fwd,
                "primer_seq_rev": rev,
                "transcript_id": g_sorted["transcript_id"].iloc[0],
                "gene_id": g_sorted["gene_id"].iloc[0],
                "gene_symbol": g_sorted["gene_symbol"].iloc[0],
                "start_up": int(g_sorted["start_up"].iloc[0]),
                "end_up": int(g_sorted["end_up"].iloc[0]),
                "start_down": int(g_sorted["start_down"].iloc[0]),
                "end_down": int(g_sorted["end_down"].iloc[0]),
                "pool_fwd": pool_fwd_vals[0],
                "pool_rev": pool_rev_vals[0],
                "num_transcripts": int(g_sorted["transcript_id"].nunique()),
                "num_genes": int(g_sorted["gene_id"].nunique()),
                "transcript_id_all": ";".join(parts),
            }
        )

    return pd.DataFrame(rows)


def _build_summary_df(
    grouped_pairs: pd.DataFrame, fwd_prefix: str, rev_prefix: str, chosen_index: int = 0
) -> pd.DataFrame:
    """Construct the final summary DataFrame matching the notebook structure.

    Adds Group/Chosen Index, geneID (without version), amplicon_index, and both the
    short primer display (L_seq/R_seq) and full order sequences including universal prefixes.
    """
    df = grouped_pairs.copy()
    if df.empty:
        return pd.DataFrame(
            columns=[
                "Group",
                "geneSymbol",
                "geneID",
                "Chosen Index",
                "amplicon_index",
                "L_seq",
                "R_seq",
                "primer_sequence_to_order_forward",
                "primer_sequence_to_order_reverse",
                "L_pool",
                "R_pool",
                "multi_mapping",
            ]
        )

    df_out = df.rename(
        columns={
            "primer_seq_fwd": "L_seq",
            "primer_seq_rev": "R_seq",
            "pool_fwd": "L_pool",
            "pool_rev": "R_pool",
            "gene_symbol": "geneSymbol",
            "transcript_id_all": "multi_mapping",
        }
    ).copy()

    df_out["Group"] = "default"
    df_out["Chosen Index"] = int(chosen_index)
    df_out["geneID"] = df_out["gene_id"].astype(str).str.split(".").str[0]

    # amplicon_index = transcript without version + _ + chosen_index
    df_out["amplicon_index"] = (
        df_out["transcript_id"].astype(str).str.split(".").str[0]
        + "_"
        + df_out["Chosen Index"].astype(str)
    )

    # full order sequences
    df_out["primer_sequence_to_order_forward"] = fwd_prefix + df_out["L_seq"].astype(
        str
    )
    df_out["primer_sequence_to_order_reverse"] = rev_prefix + df_out["R_seq"].astype(
        str
    )

    # For display L_seq/R_seq, include any tail after NNNNNN but exclude the universal prefix before NNNNNN
    fwd_tail = fwd_prefix.split("NNNNNN", 1)[-1]
    rev_tail = rev_prefix.split("NNNNNN", 1)[-1]
    df_out["L_seq"] = fwd_tail + df_out["L_seq"].astype(str)
    df_out["R_seq"] = rev_tail + df_out["R_seq"].astype(str)

    cols = [
        "Group",
        "geneSymbol",
        "geneID",
        "Chosen Index",
        "amplicon_index",
        "L_seq",
        "R_seq",
        "primer_sequence_to_order_forward",
        "primer_sequence_to_order_reverse",
        "L_pool",
        "R_pool",
        "multi_mapping",
    ]
    return df_out[cols].copy()


def recover(
    order_excel: Path | str,
    txome_fasta: Path | None = None,
    species: Literal["mouse", "human"] | None = None,
    has_linker: bool = False,
    output_order_excel: Path | str | None = None,
    output_report: Path | str | None = None,
    output_i5: Path | str | None = None,
    output_i7: Path | str | None = None,
    *,
    amplicon_length_range: tuple[int, int] = (320, 380),
    threads: int = 1,
) -> pd.DataFrame:
    """Recover the primer summary from an IDT order Excel and a transcriptome FASTA.

    This function replicates the logic in the provided notebook:
    - Parse the IDT sheet and classify primers using the known universal prefixes.
    - Locate gene-specific primer sequences across the transcriptome with `seqkit locate`.
    - Pair forward/reverse hits into amplicons, filter by length, and summarize per primer pair.
    - Optionally write the recovered summary to Excel and a multi-page PDF report.

    Args:
        txome_fasta: Path to the transcript FASTA (gz ok). If provided, overrides `species`.
        species: One of {"mouse", "human"} to resolve a packaged FASTA via `_data_dir_for_species`.
        order_excel: Path to the IDT order Excel file. First two columns must be pool and sequence.
        has_linker: Whether the panel used the extra 6bp linker after NNNNNN (ACATCA/ATAGTT).
        output_order_excel: If provided, write the recovered summary to this Excel path.
        output_report: If provided, write a PDF report with sanity-check plots.
        amplicon_length_range: Inclusive target amplicon length range used for filtering.
        threads: Number of threads for `seqkit locate` (defaults to `os.cpu_count()`).

    Returns:
        A pandas DataFrame of the recovered summary with the following columns:
        [Group, geneSymbol, geneID, Chosen Index, amplicon_index, L_seq, R_seq,
         primer_sequence_to_order_forward, primer_sequence_to_order_reverse,
         L_pool, R_pool, multi_mapping].

    Raises:
        ValueError: On invalid inputs or unexpected IDT file contents.
        FileNotFoundError: If the resolved FASTA path does not exist.
        RuntimeError: If `seqkit` is missing or fails, or if internal joins fail.
    """
    order_excel = Path(order_excel)

    # Resolve reference path
    ref_path = resolve_reference_path(txome_fasta, species)

    # Prefixes
    fwd_prefix, rev_prefix = get_prefixes(has_linker)

    # Read and classify IDT excel
    primer_info = _classify_primers(
        _read_idt_excel(order_excel), fwd_prefix, rev_prefix
    )

    # Locate patterns in the transcriptome
    primer_seqs = primer_info["primer_seq"].astype(str).str.upper().tolist()
    # Sanity check: Are primer sequences unique?
    vc = pd.Series(primer_seqs).value_counts()
    logger.info(vc[vc > 1])
    locate_df = _run_seqkit_locate(primer_seqs, ref_path, threads=threads)
    vc = pd.Series(locate_df["pattern"]).value_counts()
    logger.info(vc[vc > 1])

    # Enrich locate dataframe
    locate_df_final = _enrich_locate_df(locate_df, primer_info)
    # Sanity check: All primer sequences match their lengths.
    assert (
        (locate_df_final["end"] - locate_df_final["start"] + 1)
        == locate_df_final["primer_seq"].str.len()
    ).all()

    # Enumerate and filter amplicons
    amplicon_all = _build_amplicons(locate_df_final)
    lo, hi = amplicon_length_range
    amplicon_sub = amplicon_all[
        amplicon_all["amplicon_length"].between(lo, hi, inclusive="both")
    ].copy()

    # Group per primer pair and construct final sheet
    grouped = _group_pairs(amplicon_sub)
    # Sanity check: Ideally, it should be a one-to-one relation between forward and reverse primers.
    logger.info(
        grouped["primer_seq_fwd"].value_counts().max(),
        grouped["primer_seq_rev"].value_counts().max(),
    )
    summary_df = _build_summary_df(grouped, fwd_prefix, rev_prefix, chosen_index=0)

    # Output artifacts
    if output_order_excel is not None:
        out_xlsx = Path(output_order_excel)
        out_xlsx.parent.mkdir(parents=True, exist_ok=True)
        summary_df.to_excel(out_xlsx, index=False)

    if output_report is not None:
        make_report(
            Path(output_report),
            amplicon_all,
            amplicon_sub,
            locate_df_final,
            amplicon_length_range,
        )

    # i5/i7 short FASTAs
    _write_short_fastas(summary_df, output_i5, output_i7)

    return summary_df


def _write_short_fastas(
    summary_df: pd.DataFrame, output_i5: str | Path | None, output_i7: str | Path | None
) -> None:
    """Write i5/i7 *short* FASTAs from the recovered summary.

    The short sequences correspond to the display `L_seq`/`R_seq` columns in the
    summary (i.e., the tail after NNNNNN plus the gene-specific portion),
    mirroring the notebook behavior. FASTA record IDs are constructed from the
    unique set of gene symbols in the `multi_mapping` column for each unique
    `L_seq`/`R_seq`, joined by `|`. If an ID would collide, suffix with `.1`, `.2`, ...

    Args:
        summary_df: Final summary DataFrame from :func:`recover`.
        output_i5: Destination path for the i5 short FASTA (forward). If None, skip.
        output_i7: Destination path for the i7 short FASTA (reverse). If None, skip.
    """
    if summary_df.empty:
        return

    def _build_mapping(df: pd.DataFrame, seq_col: str) -> pd.DataFrame:
        grouped = (
            df[[seq_col, "multi_mapping"]]
            .groupby(seq_col, as_index=False)
            .agg({"multi_mapping": lambda x: ";".join(x)})
        )
        return grouped

    def _id_from_multi_mapping(multi: str) -> str:
        if not isinstance(multi, str) or not multi:
            return "NA"
        symbols: set[str] = set()
        for part in multi.split(";"):
            fields = part.split("|")
            if len(fields) >= 3 and fields[2]:
                symbols.add(fields[2])
        return "|".join(sorted(symbols)) if symbols else "NA"

    if output_i5 is not None:
        out = Path(output_i5)
        out.parent.mkdir(parents=True, exist_ok=True)
        fwd_df = _build_mapping(summary_df, "L_seq")
        used_ids: set[str] = set()
        with open(out, "w") as handle:
            for _, row in fwd_df.iterrows():
                seq: str = str(row["L_seq"]).upper()
                rid = _id_from_multi_mapping(str(row["multi_mapping"]))
                base = rid if rid else "NA"
                rid_final = base
                k = 1
                while rid_final in used_ids:
                    rid_final = f"{base}.{k}"
                    k += 1
                used_ids.add(rid_final)
                SeqIO.write(
                    SeqRecord(Seq(seq), id=rid_final, description=""), handle, "fasta"
                )

    if output_i7 is not None:
        out = Path(output_i7)
        out.parent.mkdir(parents=True, exist_ok=True)
        rev_df = _build_mapping(summary_df, "R_seq")
        used_ids: set[str] = set()
        with open(out, "w") as handle:
            for _, row in rev_df.iterrows():
                seq: str = str(row["R_seq"]).upper()
                rid = _id_from_multi_mapping(str(row["multi_mapping"]))
                base = rid if rid else "NA"
                rid_final = base
                k = 1
                while rid_final in used_ids:
                    rid_final = f"{base}.{k}"
                    k += 1
                used_ids.add(rid_final)
                SeqIO.write(
                    SeqRecord(Seq(seq), id=rid_final, description=""), handle, "fasta"
                )
