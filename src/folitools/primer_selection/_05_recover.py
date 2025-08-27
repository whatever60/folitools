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
from .recover_utils import simplify_gene_list


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setLevel(logging.INFO)  # Need both logger level and handler level to work
logger.propagate = False
logger.addHandler(handler)

SANITY_CHECK_TOP_N = 5


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
    primer_df["primer_seq_full"] = primer_df["primer_seq_full"].astype(str).str.upper()
    primer_df["primer_seq"] = primer_df["primer_seq"].astype(str).str.upper()

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

    # logger.info("Running: %s", " ".join(args))
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
    if locate_df.empty:
        # Single source of truth (local to this branch) for columns & dtypes
        spec = [
            ("primer_seq", "string"),
            ("primer_seq_full", "string"),
            ("primer_type", "string"),
            ("transcript_id", "string"),
            ("gene_id", "string"),
            ("gene_symbol", "string"),
            ("start", "Int64"),
            ("end", "Int64"),
            ("strand", "string"),
            ("pool", "string"),
        ]
        return pd.DataFrame({c: pd.Series(dtype=dt) for c, dt in spec})

    df = locate_df.copy()
    df["primer_seq"] = df["pattern"].str.upper()

    # Parse pipe-separated identifiers safely
    parts = df["seqID"].astype(str).str.split("|", expand=True)
    df["transcript_id"] = parts[0]
    df["gene_id"] = parts[1] if parts.shape[1] > 1 else "_NA"
    df["gene_symbol"] = parts[5] if parts.shape[1] > 5 else df["transcript_id"]

    merged = pd.merge(df, primer_info, how="left", on="primer_seq")
    if merged["pool"].isna().any():
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

    return locate_df_final


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
        spec = [
            ("primer_seq_fwd", "string"),
            ("primer_seq_rev", "string"),
            ("transcript_id", "string"),
            ("gene_id", "string"),
            ("gene_symbol", "string"),
            ("start_up", "Int64"),
            ("end_up", "Int64"),
            ("start_down", "Int64"),
            ("end_down", "Int64"),
            ("pool_fwd", "string"),
            ("pool_rev", "string"),
            ("flipped", "boolean"),
            ("amplicon_length", "Int64"),
        ]
        return pd.DataFrame({c: pd.Series(dtype=dt) for c, dt in spec})
    return amplicon_all


def _group_pairs(amplicon_sub: pd.DataFrame) -> pd.DataFrame:
    """Group by primer pair and summarize mappings across transcripts/genes.

    Returns a DataFrame with one row per (forward, reverse) pair with stable fields
    and a semicolon-joined "multi_mapping" column.
    """
    if amplicon_sub.empty:
        spec = [
            ("primer_seq_fwd", "string"),
            ("primer_seq_rev", "string"),
            ("transcript_id", "string"),
            ("gene_id", "string"),
            ("gene_symbol", "string"),
            ("start_up", "Int64"),
            ("end_up", "Int64"),
            ("start_down", "Int64"),
            ("end_down", "Int64"),
            ("pool_fwd", "string"),
            ("pool_rev", "string"),
            ("num_transcripts", "Int64"),
            ("num_genes", "Int64"),
            ("transcript_id_all", "string"),
        ]
        return pd.DataFrame({c: pd.Series(dtype=dt) for c, dt in spec})

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
        spec = [
            ("Group", "string"),
            ("geneSymbol", "string"),
            ("geneID", "string"),
            ("Chosen Index", "Int64"),
            ("amplicon_index", "string"),
            ("L_seq", "string"),
            ("R_seq", "string"),
            ("primer_sequence_to_order_forward", "string"),
            ("primer_sequence_to_order_reverse", "string"),
            ("L_pool", "string"),
            ("R_pool", "string"),
            ("multi_mapping", "string"),
        ]
        return pd.DataFrame({c: pd.Series(dtype=dt) for c, dt in spec})

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


def _sanity_check_primer_duplicates(primer_seqs: list[str]) -> None:
    """Sanity check: Are primer sequences unique?"""
    vc = pd.Series(primer_seqs).value_counts()
    duplicated_primers = vc[vc > 1]
    if not duplicated_primers.empty:
        logger.warning(
            "⚠️ %d duplicate primer sequences found:", len(duplicated_primers)
        )
        for seq, count in duplicated_primers.head(SANITY_CHECK_TOP_N).items():
            logger.warning("  - %s: %d occurrences", seq, count)
        if len(duplicated_primers) > SANITY_CHECK_TOP_N:
            logger.warning(
                "  ... and %d more duplicate sequences.",
                len(duplicated_primers) - SANITY_CHECK_TOP_N,
            )
    else:
        logger.info("✅ All primer sequences are unique.")


def _sanity_check_primer_duplicates_with_type(primer_info: pd.DataFrame) -> None:
    """Sanity check: Are primer sequences unique with respect to primer type?"""
    duplicate_counts = (
        primer_info.groupby(["primer_type", "primer_seq"])
        .size()
        .reset_index(name="counts")
    )
    duplicate_counts = duplicate_counts[duplicate_counts["counts"] > 1]

    if duplicate_counts.empty:
        logger.info("✅ All primer sequences are unique with respect to primer type.")
    else:
        logger.warning(
            "⚠️ %d duplicate primer sequences with same type found:",
            len(duplicate_counts),
        )
        for _, row in duplicate_counts.head(SANITY_CHECK_TOP_N).iterrows():
            logger.warning(
                "  - (%s, %s): %d occurrences",
                row["primer_type"],
                row["primer_seq"],
                row["counts"],
            )
        if len(duplicate_counts) > SANITY_CHECK_TOP_N:
            logger.warning(
                "  ... and %d more duplicate sequences with same type.",
                len(duplicate_counts) - 10,
            )


def _sanity_check_seqkit_patterns(locate_df: pd.DataFrame) -> None:
    """Sanity check: Number of matches in the transcriptome for each primer. We could
    see some primers have thousands of matches.
    """
    vc = locate_df["pattern"].value_counts()
    if vc.empty:
        logger.info("⚪ No seqkit patterns to check.")
        return

    multi_match_patterns = vc[vc > 1]
    if multi_match_patterns.empty:
        logger.info("✅ All primers have one single match in the transcriptome.")
    else:
        logger.warning(
            "⚠️ %d primers have multiple matches in the transcriptome:",
            len(multi_match_patterns),
        )
        for pattern, count in multi_match_patterns.head(SANITY_CHECK_TOP_N).items():
            logger.warning("  - %s: %d matches", pattern, count)
        if len(multi_match_patterns) > SANITY_CHECK_TOP_N:
            logger.warning(
                "  ... and %d more patterns with multiple matches.",
                len(multi_match_patterns) - SANITY_CHECK_TOP_N,
            )


def _sanity_check_multi_location_binding(locate_df_final: pd.DataFrame) -> None:
    """Sanity check: Ideally, no duplicate primer_seq and transcript_id pairs. In reality,
    a primer sequence can bind to multiple locations on the same transcript if the
    transcript sequence is repetitive.
    """
    binding_counts = locate_df_final.groupby(["primer_seq", "transcript_id"]).size()
    multi_bindings = binding_counts[binding_counts > 1].sort_values(ascending=False)

    if multi_bindings.empty:
        logger.info("✅ No primers bind to multiple locations on the same transcript.")
    else:
        logger.warning(
            "⚠️ Found %d instances of primers binding to multiple locations on the same transcript:",
            len(multi_bindings),
        )
        for (
            (primer_seq, transcript_id),
            count,
        ) in multi_bindings.head(SANITY_CHECK_TOP_N).items():
            logger.warning(
                "  - Primer %s binds to %d locations on transcript %s",
                primer_seq,
                count,
                transcript_id,
            )
        if len(multi_bindings) > SANITY_CHECK_TOP_N:
            logger.warning(
                "  ... and %d more multi-location binding instances.",
                len(multi_bindings) - SANITY_CHECK_TOP_N,
            )


def _sanity_check_multi_gene_binding(locate_df_final: pd.DataFrame) -> None:
    """Sanity check: Ideally, a primer_seq matches to one gene_id only. In reality, a
    primer sequence could bind to multiple locations on multiple transcripts of the same gene.
    """
    gene_counts = locate_df_final.groupby("primer_seq")["gene_id"].nunique()
    multi_gene_bindings = gene_counts[gene_counts > 1].sort_values(ascending=False)

    if multi_gene_bindings.empty:
        logger.info("✅ All primers bind to a single gene.")
    else:
        logger.warning(
            "⚠️ Found %d primers binding to multiple genes:", len(multi_gene_bindings)
        )
        for primer_seq, num_genes in multi_gene_bindings.head(
            SANITY_CHECK_TOP_N
        ).items():
            genes = locate_df_final[locate_df_final["primer_seq"] == primer_seq][
                "gene_id"
            ].unique()
            # genes_str = ", ".join(genes[:5])
            if len(genes) > SANITY_CHECK_TOP_N:
                genes = genes[:SANITY_CHECK_TOP_N].tolist() + [
                    f"... and {len(genes) - SANITY_CHECK_TOP_N} more genes"
                ]
            else:
                genes = genes.tolist()
            logger.warning(
                "  - Primer %s binds to %d genes: %s",
                primer_seq,
                num_genes,
                ", ".join(genes),
            )
        if len(multi_gene_bindings) > SANITY_CHECK_TOP_N:
            logger.warning(
                "  ... and %d more multi-gene binding primers.",
                len(multi_gene_bindings) - SANITY_CHECK_TOP_N,
            )


def _sanity_check_sequence_lengths(locate_df_final: pd.DataFrame) -> None:
    """Sanity check: All primer sequences match their lengths."""
    mismatches = locate_df_final[
        (locate_df_final["end"] - locate_df_final["start"] + 1)
        != locate_df_final["primer_seq"].str.len()
    ]
    if mismatches.empty:
        logger.info("✅ All primer sequence lengths match their locate positions.")
    else:
        logger.error(
            f"❌ {len(mismatches)} primers have length mismatches in locate results!"
        )
        raise ValueError(
            f"{len(mismatches)} primers have length mismatches in locate results!"
        )


def _sanity_check_genes_per_primer_pair(amplicon_sub: pd.DataFrame) -> None:
    """Sanity check: Check the number of genes amplified by each primer pair."""
    if amplicon_sub.empty:
        logger.info("⚪ No amplicons to analyze for gene specificity.")
        return

    genes_per_pair = amplicon_sub.groupby(["primer_seq_fwd", "primer_seq_rev"])[
        "gene_id"
    ].nunique()
    multi_gene_pairs = genes_per_pair[genes_per_pair > 1].sort_values(ascending=False)

    if multi_gene_pairs.empty:
        logger.info("✅ All primer pairs amplify a single gene.")
    else:
        logger.warning(
            "⚠️ %d primer pairs amplify more than one gene:",
            len(multi_gene_pairs),
        )
        for (fwd, rev), count in multi_gene_pairs.head(SANITY_CHECK_TOP_N).items():
            logger.warning("  - %s/%s: %d genes", fwd, rev, count)
        if len(multi_gene_pairs) > SANITY_CHECK_TOP_N:
            logger.warning(
                "  ... and %d more multi-gene pairs.",
                len(multi_gene_pairs) - SANITY_CHECK_TOP_N,
            )


def _sanity_check_amplicons_per_pair(amplicon_sub: pd.DataFrame) -> None:
    """Sanity check: Number of amplicons formed by each primer pair."""
    if amplicon_sub.empty:
        logger.info("⚪ No amplicons to analyze.")
        return

    amplicons_per_pair = amplicon_sub.groupby(
        ["primer_seq_fwd", "primer_seq_rev"]
    ).size()
    multi_amplicons = amplicons_per_pair[amplicons_per_pair > 1].sort_values(
        ascending=False
    )

    if multi_amplicons.empty:
        logger.info("✅ All primer pairs form a single amplicon.")
    else:
        logger.warning(
            "⚠️ %d primer pairs form multiple amplicons:",
            len(multi_amplicons),
        )
        for (fwd, rev), count in multi_amplicons.head(SANITY_CHECK_TOP_N).items():
            logger.warning("  - %s/%s: %d amplicons", fwd, rev, count)
        if len(multi_amplicons) > SANITY_CHECK_TOP_N:
            logger.warning(
                "  ... and %d more multi-amplicon pairs.",
                len(multi_amplicons) - SANITY_CHECK_TOP_N,
            )


def _sanity_check_cross_pool_amplicons(amplicon_sub: pd.DataFrame) -> None:
    """Sanity check: Cross-pool amplicon (i.e., the two primers of the amplicon are in
    different pools)
    """
    if amplicon_sub.empty:
        logger.info("⚪ No amplicons to check for cross-pool formation.")
        return

    cross_pool_mask = amplicon_sub["pool_fwd"] != amplicon_sub["pool_rev"]
    cross_pool_count = cross_pool_mask.sum()

    if cross_pool_count == 0:
        logger.info("✅ No cross-pool amplicons were formed.")
    else:
        logger.warning(
            "⚠️ Found %d cross-pool amplicons (%.1f%% of total).",
            cross_pool_count,
            (cross_pool_count / len(amplicon_sub)) * 100,
        )


def _sanity_check_primer_pair_relationships(grouped: pd.DataFrame) -> None:
    """Sanity check: Ideally, it should be a one-to-one relation between forward and
    reverse primers.
    """
    if grouped.empty:
        logger.info("⚪ No primer pairs to analyze for relationships.")
        return

    fwd_counts = grouped["primer_seq_fwd"].value_counts()
    rev_counts = grouped["primer_seq_rev"].value_counts()
    fwd_max = fwd_counts.max()
    rev_max = rev_counts.max()

    if fwd_max <= 1 and rev_max <= 1:
        logger.info("✅ One-to-one primer pair relationships maintained.")
    else:
        logger.warning("⚠️ Some primers are reused in multiple pairs:")
        if fwd_max > 1:
            logger.warning("  - Forward primers are used up to %d times.", fwd_max)
        if rev_max > 1:
            logger.warning("  - Reverse primers are used up to %d times.", rev_max)


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
    simplify_gene_name: bool = True,
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
            Use -1 to indicate no bound (e.g., (-1, 380) for no minimum, (320, -1) for no maximum).
        threads: Number of threads for `seqkit locate` (defaults to `os.cpu_count()`).
        simplify_gene_name: Whether to simplify gene names by removing redundant/uninformative tokens (default: True).

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
    _sanity_check_primer_duplicates(primer_info["primer_seq"].tolist())
    _sanity_check_primer_duplicates_with_type(primer_info)

    primer_info = primer_info.drop_duplicates(["primer_seq", "primer_type"])
    primer_seqs_fwd = primer_info.query("primer_type == 'fwd'")["primer_seq"].tolist()
    primer_seqs_rev = primer_info.query("primer_type == 'rev'")["primer_seq"].tolist()

    locate_df = _run_seqkit_locate(
        primer_seqs_fwd + primer_seqs_rev, ref_path, threads=threads
    )
    _sanity_check_seqkit_patterns(locate_df)

    # Enrich locate dataframe
    locate_df_final = _enrich_locate_df(locate_df, primer_info)
    _sanity_check_multi_location_binding(locate_df_final)
    _sanity_check_multi_gene_binding(locate_df_final)
    _sanity_check_sequence_lengths(locate_df_final)

    # Enumerate and filter amplicons
    amplicon_all = _build_amplicons(locate_df_final)

    lo, hi = amplicon_length_range
    
    # Handle -1 as infinite bounds for filtering
    if lo == -1 and hi == -1:
        # No filtering if both bounds are -1
        amplicon_sub = amplicon_all.copy()
    elif lo == -1:
        # No lower bound, only upper bound
        amplicon_sub = amplicon_all[amplicon_all["amplicon_length"] <= hi].copy()
    elif hi == -1:
        # No upper bound, only lower bound
        amplicon_sub = amplicon_all[amplicon_all["amplicon_length"] >= lo].copy()
    else:
        # Both bounds specified
        amplicon_sub = amplicon_all[
            amplicon_all["amplicon_length"].between(lo, hi, inclusive="both")
        ].copy()
    _sanity_check_genes_per_primer_pair(amplicon_sub)
    _sanity_check_amplicons_per_pair(amplicon_sub)
    _sanity_check_cross_pool_amplicons(amplicon_sub)
    if output_report is not None:
        make_report(
            Path(output_report),
            amplicon_all,
            amplicon_sub,
            locate_df_final,
            amplicon_length_range,
        )

    # Group per primer pair and construct final sheet
    grouped = _group_pairs(amplicon_sub)
    _sanity_check_primer_pair_relationships(grouped)

    summary_df = _build_summary_df(grouped, fwd_prefix, rev_prefix, chosen_index=0)

    # Output artifacts
    if output_order_excel is not None:
        out_xlsx = Path(output_order_excel)
        out_xlsx.parent.mkdir(parents=True, exist_ok=True)
        summary_df.to_excel(out_xlsx, index=False)

    # Create SeqRecord lists for i5/i7 short FASTAs
    i5_records = _create_fasta_records(
        summary_df.set_index("L_seq")["multi_mapping"],
        primer_seqs_fwd,
        simplify_gene_name=simplify_gene_name,
    )
    i7_records = _create_fasta_records(
        summary_df.set_index("R_seq")["multi_mapping"],
        primer_seqs_rev,
        simplify_gene_name=simplify_gene_name,
    )

    # Write i5/i7 short FASTAs if output paths are provided
    if output_i5 is not None:
        SeqIO.write(i5_records, str(output_i5), "fasta")
    if output_i7 is not None:
        SeqIO.write(i7_records, str(output_i7), "fasta")

    return summary_df


def collapse_to_id(series: pd.Series) -> str:
    symbols: set[str] = set()
    # Gather unique symbols from "multi_mapping" entries across the group
    for multi in series.dropna().astype(str):
        for part in multi.split(";"):
            fields = part.split("|")
            assert len(fields) >= 3
            symbols.add(fields[2])
    rid = "|".join(sorted(symbols)) if symbols else "_NA"
    # Keep existing behavior of simplifying the gene list
    return rid


def _create_fasta_records(
    primer_mappings: pd.Series, seqs_all: list[str], *, simplify_gene_name: bool = True
) -> list[SeqRecord]:
    """Create a list of SeqRecord objects from the summary DataFrame.

    The short sequences correspond to the display `L_seq`/`R_seq` columns in the
    summary. FASTA record IDs are constructed from the unique set of gene symbols
    in the `multi_mapping` column for each unique sequence, joined by `|`.
    If an ID would collide, a numeric suffix is added.

    Args:
        primer_mappings: A pandas Series with sequences as index (can have duplicates)
                         and "multi_mapping" strings as values.
        seqs_all: A list of all expected sequences ("L_seq" or "R_seq").
        simplify_gene_name: Whether to simplify gene names by removing redundant/uninformative tokens.

    Returns:
        A list of SeqRecord objects.
    """
    # Add missing sequences with _NA mapping
    present_seqs = set(primer_mappings.index)
    missing_seqs = [s for s in seqs_all if s not in present_seqs]
    if missing_seqs:
        missing_series = pd.Series("_NA|_NA|_NA", index=missing_seqs)
        primer_mappings = pd.concat([primer_mappings, missing_series])

    # Collapse multi-mapping to a simplified record_id per unique sequence.
    # Group by index (sequence) and apply collapse_to_id
    mapping_df = (
        primer_mappings.groupby(primer_mappings.index)
        .apply(
            lambda series: "|".join(
                simplify_gene_list(
                    collapse_to_id(series).split("|"), 
                    collapse_families=True
                )
            )
            if simplify_gene_name
            else collapse_to_id(series)
        )
        .reset_index()
    )
    mapping_df.columns = ["seq", "record_id"]

    used_ids: set[str] = set()
    records: list[SeqRecord] = []

    for _, row in mapping_df.iterrows():
        seq: str = str(row["seq"]).upper()
        rid: str = str(row["record_id"])
        rid_final = rid
        k = 1
        while rid_final in used_ids:
            rid_final = f"{rid}.{k}"
            k += 1
        used_ids.add(rid_final)
        records.append(SeqRecord(Seq(seq), id=rid_final, description=""))

    return records
