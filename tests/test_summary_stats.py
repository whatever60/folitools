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


def _write_add_tags_log(path: Path, rows: list[tuple]) -> None:
    """Rows are ``(sample, total_r1, not_na_adapter, good_umi)`` and may
    optionally append ``(mapped, assigned)`` and/or ``(counted,
    counted_assigned)``. Trailing fields must come in pairs in that
    order, matching the SUMMARY line written by ``foli_add_tags --log``.
    """
    lines: list[str] = []
    for row in rows:
        if len(row) not in (4, 6, 8):
            raise AssertionError(f"unexpected SUMMARY row arity: {row}")
        s, t, n, g, *rest = row
        parts = [
            f"cell_tag={s}",
            f"total_r1={t}",
            f"not_na_adapter={n}",
            f"good_umi={g}",
        ]
        if len(rest) >= 2:
            parts += [f"mapped={rest[0]}", f"assigned={rest[1]}"]
        if len(rest) >= 4:
            parts += [f"counted={rest[2]}", f"counted_assigned={rest[3]}"]
        lines.append("SUMMARY " + " ".join(parts) + "\n")
    path.write_text("".join(lines))


def _write_group_tsv(path: Path, rows: list[tuple[str, str]]) -> None:
    """Rows are ``(read_id, gene)``. Other columns are stubbed."""
    import gzip

    header = (
        "read_id\tcontig\tposition\tgene\tumi\tumi_count\t"
        "final_umi\tfinal_umi_count\tunique_id\n"
    )
    body = "".join(
        f"{rid}\tchr1\t1\t{gene}\tAAAAAA\t1\tAAAAAA\t1\t{i}\n"
        for i, (rid, gene) in enumerate(rows)
    )
    payload = (header + body).encode()
    path.parent.mkdir(parents=True, exist_ok=True)
    if str(path).endswith(".gz"):
        with gzip.open(path, "wb") as fh:
            fh.write(payload)
    else:
        path.write_bytes(payload)


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

    # add_tags SUMMARY: (sample, total_r1, not_na_adapter, good_umi, mapped,
    # assigned, counted, counted_assigned). counted = mapped & good_umi;
    # counted_assigned = counted & assigned. Row sum of raw matrix must
    # equal counted_assigned.
    _write_add_tags_log(
        tmp_path / "s1.add_tags.log",
        [("s1", 900, 800, 700, 840, 500, 690, 400)],
    )
    _write_add_tags_log(
        tmp_path / "s2.add_tags.log",
        [("s2", 760, 650, 550, 710, 400, 540, 280)],
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
    # Columns: raw, qc, long, not_na, good_umi, mapped, assigned,
    #          counted, counted_assigned, n_umi, n_genes
    assert result.loc["s1"].tolist() == [1000, 950, 900, 800, 700, 840, 500, 690, 400, 280, 2]
    assert result.loc["s2"].tolist() == [800, 780, 760, 650, 550, 710, 400, 540, 280, 150, 1]


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


def test_summary_stats_rejects_dag_violation(tmp_path: Path) -> None:
    _write_seqkit_stats(
        tmp_path / "fastq.stats",
        [("s1_S1_R1_001.fastq.gz", 100), ("s1_S1_R2_001.fastq.gz", 100)],
    )
    _write_seqkit_stats(
        tmp_path / "fastp.stats",
        [("s1_1.fq.gz", 200), ("s1_2.fq.gz", 200)],
    )

    with pytest.raises(AssertionError, match="DAG violated"):
        summary_stats(
            fastq_stats=str(tmp_path / "fastq.stats"),
            fastp_stats=str(tmp_path / "fastp.stats"),
        )


def test_summary_stats_rejects_mapping_branch_violation(tmp_path: Path) -> None:
    """assigned > mapped on the bottom branch should fire the DAG check."""
    # add_tags log only — emit mapped < assigned to violate the mapped→assigned edge.
    _write_add_tags_log(
        tmp_path / "s1.add_tags.log",
        [("s1", 900, 800, 700, 410, 500)],
    )

    with pytest.raises(AssertionError, match=r"mapped_depth.*assigned_depth"):
        summary_stats(add_tags_logs=str(tmp_path / "*.add_tags.log"))


def test_summary_stats_strict_passes_when_sources_agree(tmp_path: Path) -> None:
    """strict=True succeeds when add_tags counts == TSV row counts == matrix row sum."""
    _write_add_tags_log(
        tmp_path / "s1.add_tags.log",
        [("s1", 900, 800, 700, 840, 500, 4, 3)],
    )
    _write_group_tsv(
        tmp_path / "s1.group.tsv.gz",
        [
            ("r1", "GENE1,FWD+REV"),
            ("r2", "GENE2,FWD+REV"),
            ("r3", "GENE1,FWD+REV"),
            ("r4", "Unassigned,FWD+REV"),
        ],
    )
    raw = pd.DataFrame({"GENE1": [2], "GENE2": [1]}, index=["s1"])
    raw.index.name = "sample"
    _write_count_matrix(tmp_path / "raw.tsv", raw)

    result = summary_stats(
        add_tags_logs=str(tmp_path / "*.add_tags.log"),
        group_tsvs=str(tmp_path / "*.group.tsv.gz"),
        count_matrix_raw=str(tmp_path / "raw.tsv"),
        strict=True,
    )

    assert result.loc["s1", "counted_depth"] == 4
    assert result.loc["s1", "counted_assigned_depth"] == 3


def test_summary_stats_strict_rejects_count_disagreement(tmp_path: Path) -> None:
    """strict=True asserts when add_tags counted != group.tsv.gz row count."""
    _write_add_tags_log(
        tmp_path / "s1.add_tags.log",
        [("s1", 900, 800, 700, 840, 500, 5, 3)],
    )
    _write_group_tsv(
        tmp_path / "s1.group.tsv.gz",
        [("r1", "GENE1,FWD+REV"), ("r2", "GENE2,FWD+REV")],
    )

    with pytest.raises(AssertionError, match="counted_depth=5.*group.tsv.gz row count"):
        summary_stats(
            add_tags_logs=str(tmp_path / "*.add_tags.log"),
            group_tsvs=str(tmp_path / "*.group.tsv.gz"),
            strict=True,
        )


def test_summary_stats_strict_rejects_matrix_disagreement(tmp_path: Path) -> None:
    """strict=True asserts when counted_assigned != row sum of raw matrix."""
    _write_add_tags_log(
        tmp_path / "s1.add_tags.log",
        [("s1", 900, 800, 700, 840, 500, 690, 400)],
    )
    raw = pd.DataFrame({"GENE1": [200], "GENE2": [199]}, index=["s1"])
    raw.index.name = "sample"
    _write_count_matrix(tmp_path / "raw.tsv", raw)

    with pytest.raises(AssertionError, match=r"counted_assigned_depth=400.*count_matrix_raw"):
        summary_stats(
            add_tags_logs=str(tmp_path / "*.add_tags.log"),
            count_matrix_raw=str(tmp_path / "raw.tsv"),
            strict=True,
        )


def test_summary_stats_rejects_unnamed_add_tags_log(tmp_path: Path) -> None:
    _write_add_tags_log(
        tmp_path / "unknown.add_tags.log",
        [("-", 100, 80, 70)],
    )
    with pytest.raises(ValueError, match="cell_tag=-"):
        summary_stats(add_tags_logs=str(tmp_path / "*.add_tags.log"))


def test_summary_cli_writes_table(tmp_path: Path) -> None:
    from folitools.cli import summary as summary_cmd

    _write_seqkit_stats(
        tmp_path / "fastq.stats",
        [("s1_S1_R1_001.fastq.gz", 100), ("s1_S1_R2_001.fastq.gz", 100)],
    )
    _write_seqkit_stats(
        tmp_path / "fastp.stats",
        [("s1_1.fq.gz", 80), ("s1_2.fq.gz", 80)],
    )

    output = tmp_path / "summary.tsv"
    summary_cmd(
        output=str(output),
        fastq_stats=str(tmp_path / "fastq.stats"),
        fastp_stats=str(tmp_path / "fastp.stats"),
    )

    df = pd.read_csv(output, sep="\t", index_col=0)
    assert df.index.name.startswith("folitools ")
    assert list(df.columns) == list(METRIC_COLUMNS)
    assert df.loc["s1", "raw_depth"] == 100
    assert df.loc["s1", "pass_qc_depth"] == 80
