"""Tests for folitools.summary.summary_stats."""

from pathlib import Path
import textwrap

import pandas as pd
import pytest

from folitools.summary import METRIC_COLUMNS, summary_stats


def _write_seqkit_stats(path: Path, rows: list[tuple[str, int]]) -> None:
    header = (
        "file\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len\t"
        "Q1\tQ2\tQ3\tsum_gap\tN50\tN50_num\tQ20(%)\tQ30(%)\tAvgQual\tGC(%)\tsum_n\n"
    )
    path.write_text(
        header
        + "".join(
            f"{fname}\tFASTQ\tDNA\t{n}\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n"
            for fname, n in rows
        )
    )


def _write_star_log(path: Path, input_reads: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        textwrap.dedent(
            f"""\
                                     Started job on |	Dec 02 06:23:12
                              Number of input reads |	{input_reads}
                          Average input read length |	288
            """
        )
    )


def _write_add_tags_log(
    path: Path, rows: list[tuple[str, int, int, int]]
) -> None:
    path.write_text(
        "".join(
            f"SUMMARY cell_tag={s} total_r1={t} good_umi={g} not_na_adapter={a}\n"
            for s, t, g, a in rows
        )
    )


def _write_count_matrix(path: Path, df: pd.DataFrame) -> None:
    df.to_csv(path, sep="\t", index_label="folitools 0.0.0")


def test_summary_stats_end_to_end(tmp_path: Path) -> None:
    _write_seqkit_stats(
        tmp_path / "fastq.stats",
        [
            ("s1_S1_R1_001.fastq.gz", 1000),
            ("s1_S1_R2_001.fastq.gz", 1000),
            ("s2_S2_R1_001.fastq.gz", 800),
            ("s2_S2_R2_001.fastq.gz", 800),
        ],
    )
    _write_seqkit_stats(
        tmp_path / "fastp.stats",
        [
            ("s1_1.fq.gz", 950),
            ("s1_2.fq.gz", 950),
            ("s2_1.fq.gz", 780),
            ("s2_2.fq.gz", 780),
        ],
    )
    _write_star_log(tmp_path / "star" / "s1" / "Log.final.out", 900)
    _write_star_log(tmp_path / "star" / "s2" / "Log.final.out", 760)

    _write_add_tags_log(
        tmp_path / "s1.add_tags.log",
        [("s1", 900, 800, 700)],
    )
    _write_add_tags_log(
        tmp_path / "s2.add_tags.log",
        [("s2", 760, 650, 550)],
    )

    raw = pd.DataFrame(
        {"GENE1": [300, 200], "GENE2": [100, 80]}, index=["s1", "s2"]
    )
    raw.index.name = "sample"
    _write_count_matrix(tmp_path / "raw.tsv", raw)

    dedup = pd.DataFrame(
        {"GENE1": [200, 150], "GENE2": [80, 0]}, index=["s1", "s2"]
    )
    dedup.index.name = "sample"
    _write_count_matrix(tmp_path / "dedup.tsv", dedup)

    result = summary_stats(
        fastq_stats=str(tmp_path / "fastq.stats"),
        fastp_stats=str(tmp_path / "fastp.stats"),
        star_logs=str(tmp_path / "star" / "*" / "Log.final.out"),
        add_tags_logs=str(tmp_path / "*.add_tags.log"),
        count_matrix_raw=str(tmp_path / "raw.tsv"),
        count_matrix_dedup=str(tmp_path / "dedup.tsv"),
    )

    assert list(result.columns) == list(METRIC_COLUMNS)
    assert list(result.index) == ["s1", "s2"]
    assert result.loc["s1"].tolist() == [1000, 950, 900, 800, 700, 400, 280, 2]
    assert result.loc["s2"].tolist() == [800, 780, 760, 650, 550, 280, 150, 1]


def test_summary_stats_missing_sources_yield_nan(tmp_path: Path) -> None:
    _write_seqkit_stats(
        tmp_path / "fastq.stats",
        [("s1_S1_R1_001.fastq.gz", 100), ("s1_S1_R2_001.fastq.gz", 100)],
    )

    result = summary_stats(fastq_stats=str(tmp_path / "fastq.stats"))

    assert result.loc["s1", "raw_depth"] == 100
    for col in METRIC_COLUMNS:
        if col == "raw_depth":
            continue
        assert pd.isna(result.loc["s1", col]), f"{col} should be NA"


def test_summary_stats_rejects_non_monotonic(tmp_path: Path) -> None:
    _write_seqkit_stats(
        tmp_path / "fastq.stats",
        [("s1_S1_R1_001.fastq.gz", 100), ("s1_S1_R2_001.fastq.gz", 100)],
    )
    _write_seqkit_stats(
        tmp_path / "fastp.stats",
        [("s1_1.fq.gz", 200), ("s1_2.fq.gz", 200)],
    )

    with pytest.raises(AssertionError, match="monotonically non-increasing"):
        summary_stats(
            fastq_stats=str(tmp_path / "fastq.stats"),
            fastp_stats=str(tmp_path / "fastp.stats"),
        )


def test_summary_stats_rejects_unnamed_add_tags_log(tmp_path: Path) -> None:
    _write_add_tags_log(
        tmp_path / "unknown.add_tags.log",
        [("-", 100, 80, 70)],
    )
    with pytest.raises(ValueError, match="cell_tag=-"):
        summary_stats(add_tags_logs=str(tmp_path / "*.add_tags.log"))
