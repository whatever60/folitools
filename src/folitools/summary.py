"""Per-sample pipeline-stage counts for foli runs.

``summary_stats`` aggregates read counts from each stage of the foli pipeline
into a single DataFrame (rows = samples, columns = metrics). Any source left
as ``None`` yields an all-NaN column.

The metrics form a DAG, not a single chain: from the post-QC long-read pool
two parallel funnels merge at ``counted_depth`` / ``counted_assigned_depth``.

* Library-quality funnel: ``long_read_depth`` → ``not_na_adapter_depth``
  (primer pair recognized) → ``good_umi_depth`` (also a clean UMI; the
  reads that carry a ``UC`` tag).
* Mapping/annotation funnel: ``long_read_depth`` → ``mapped_depth``
  (both primary mates aligned) → ``assigned_depth`` (at least one mate
  carries a real gene id in ``XT``).

The two funnels then merge:

* ``counted_depth`` = ``mapped`` ∩ ``good_umi``: the QNAMEs that umi_tools
  emits to ``group.tsv.gz`` under foli's current options
  (``--chimeric-pairs use``, ``--unmapped-reads discard``,
  ``--unpaired-reads use``); chimeric pairs are kept as long as both ends
  are mapped and the ``UC`` tag is present.
* ``counted_assigned_depth`` = ``counted`` ∩ ``assigned``: also requires
  the QNAME's ``XF`` not to start with ``Unassigned,``, mirroring the
  filter inside ``folitools.get_matrix``. Equals the row sum of the
  pre-dedup count matrix.

All six ``add_tags``-derived metrics are counted per-QNAME inside
``foli_add_tags`` and emitted on its ``SUMMARY`` line, so they share
units with the rest of the table — unlike featureCounts' ``*.summary``
``Assigned`` row, which counts reads (R1+R2) under foli's ``-p`` (no
``--countReadPairs``) invocation and is therefore not used here.

Each edge of the DAG must be non-increasing (parent ≥ child) on every
present-value pair; an ``AssertionError`` is raised otherwise so pipeline
regressions surface immediately. Edges with a NaN endpoint are skipped.

For reproducibility audits: every log file this parser reads — and every
log this pipeline produces more broadly — also carries the folitools
version that wrote it. JSON logs (fastp, cutadapt) carry a top-level
``folitools_version`` key; plain-text logs (foli_add_tags, umi_tools
group) carry a trailing ``# folitools <version>`` line; parquet outputs
(``foli get-read-stats``) carry it in schema key-value metadata; count
matrices and the summary table carry it in the index label. The seqkit
``<dir>.stats`` TSV is the lone exception (decorating it would break
parsing) — pair it with the matching JSON for that step.
"""

from pathlib import Path
import re

import pandas as pd

from .utils import expand_path_to_list


METRIC_COLUMNS = (
    "raw_depth",
    "pass_qc_depth",
    "long_read_depth",
    "not_na_adapter_depth",
    "good_umi_depth",
    "mapped_depth",
    "assigned_depth",
    "counted_depth",
    "counted_assigned_depth",
    "n_umi",
    "n_genes",
)

# Edges of the DAG used by `_assert_dag_non_increasing`. Parent ≥ child must
# hold whenever both endpoints are present. The structure encodes:
#   raw → qc → long ─┬→ not_na_adapter → good_umi ─┐
#                    │                              ↓
#                    │                            counted → counted_assigned → n_umi → n_genes
#                    │                              ↑              ↑
#                    └→ mapped ─┬────────────────────┘              │
#                               └→ assigned ─────────────────────────┘
_DAG_EDGES: tuple[tuple[str, str], ...] = (
    ("raw_depth", "pass_qc_depth"),
    ("pass_qc_depth", "long_read_depth"),
    ("long_read_depth", "not_na_adapter_depth"),
    ("not_na_adapter_depth", "good_umi_depth"),
    ("long_read_depth", "mapped_depth"),
    ("mapped_depth", "assigned_depth"),
    ("good_umi_depth", "counted_depth"),
    ("mapped_depth", "counted_depth"),
    ("counted_depth", "counted_assigned_depth"),
    ("assigned_depth", "counted_assigned_depth"),
    ("counted_assigned_depth", "n_umi"),
    ("n_umi", "n_genes"),
)

_STAR_INPUT_RE = re.compile(r"Number of input reads\s*\|\s*(\d+)")
_SUMMARY_KV_RE = re.compile(r"(\w+)=(\S+)")


def _sample_from_seqkit_file(file_field: str) -> str:
    """Replicates ``quality.read_stat`` sample extraction so stats sources align."""
    return Path(file_field).name.split(".")[0].split("_")[0]


def _is_r1(file_field: str) -> bool:
    return "R1_001" in file_field or "_1." in file_field


def _read_seqkit_r1_counts(path: str) -> pd.Series:
    """Extract per-sample R1 num_seqs from a seqkit --tabular stats file."""
    df = pd.read_table(path)
    df = df[df["file"].map(_is_r1)]
    df = df.assign(sample=df["file"].map(_sample_from_seqkit_file))
    return df.set_index("sample")["num_seqs"].astype("int64")


def _read_star_input_reads(paths: list[str]) -> pd.Series:
    """Parse STAR ``Log.final.out`` files for ``Number of input reads``.

    Sample name is taken from the parent directory, matching the
    ``star_foli/<sample>/Log.final.out`` layout written by ``foli map``.
    """
    result: dict[str, int] = {}
    for p in paths:
        sample = Path(p).parent.name
        with open(p) as fh:
            for line in fh:
                if (m := _STAR_INPUT_RE.search(line)) is not None:
                    result[sample] = int(m.group(1))
                    break
    return pd.Series(result, dtype="int64")


_ADD_TAGS_FIELDS = (
    "not_na_adapter",
    "good_umi",
    "mapped",
    "assigned",
    "counted",
    "counted_assigned",
)


def _read_add_tags_summary(paths: list[str]) -> dict[str, pd.Series]:
    """Parse SUMMARY lines from ``foli_add_tags --log`` outputs.

    Returns a dict keyed by metric name. Every value is a per-QNAME
    counter incremented inside ``foli_add_tags`` so they share units
    across the row. Older SUMMARY lines that predate ``counted`` /
    ``counted_assigned`` simply omit those keys, yielding empty series.
    """
    out: dict[str, dict[str, int]] = {k: {} for k in _ADD_TAGS_FIELDS}
    for p in paths:
        with open(p) as fh:
            for line in fh:
                if not line.startswith("SUMMARY "):
                    continue
                fields = dict(_SUMMARY_KV_RE.findall(line))
                sample = fields.get("cell_tag", "-")
                if sample == "-":
                    raise ValueError(
                        f"SUMMARY line in {p} has cell_tag=-; rerun "
                        f"foli_add_tags with --cell_tag so samples are identifiable"
                    )
                for key in _ADD_TAGS_FIELDS:
                    if key in fields:
                        out[key][sample] = int(fields[key])
    return {k: pd.Series(v, dtype="int64") for k, v in out.items()}


def _read_count_matrix(path: str) -> pd.DataFrame:
    """Load a ``foli get-count-mtx`` matrix (samples x genes) into a DataFrame."""
    sep = "," if path.endswith((".csv", ".csv.gz")) else "\t"
    df = pd.read_csv(path, sep=sep, index_col=0)
    df.index.name = "sample"
    return df


def _sample_from_group_tsv(path: str) -> str:
    """``<sample>.group.tsv.gz`` → ``<sample>``."""
    name = Path(path).name
    for suffix in (".group.tsv.gz", ".group.tsv", ".tsv.gz", ".tsv"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return Path(path).stem


def _read_group_tsv_counts(paths: list[str]) -> tuple[pd.Series, pd.Series]:
    """Count rows in each ``group.tsv.gz`` (umi_tools group output).

    Returns ``(total, assigned)`` per sample where ``total`` is the
    number of data rows and ``assigned`` is the number of rows whose
    ``gene`` column does not start with ``Unassigned,`` (mirroring the
    filter in :mod:`folitools.get_matrix`).
    """
    total: dict[str, int] = {}
    assigned: dict[str, int] = {}
    for p in paths:
        sample = _sample_from_group_tsv(p)
        # polars handles .gz transparently and is fast for row counting.
        # Avoid loading the full table into memory by streaming.
        import polars as pl

        scan = pl.scan_csv(p, separator="\t")
        total[sample] = int(scan.select(pl.len()).collect()["len"][0])
        assigned[sample] = int(
            scan.filter(~pl.col("gene").str.starts_with("Unassigned,"))
            .select(pl.len())
            .collect()["len"][0]
        )
    return (
        pd.Series(total, dtype="int64"),
        pd.Series(assigned, dtype="int64"),
    )


def _as_path_list(arg: str | list[str] | None, suffix: str | None) -> list[str]:
    if arg is None:
        return []
    return expand_path_to_list(arg, suffix=suffix or "")


def _assert_dag_non_increasing(df: pd.DataFrame) -> None:
    """Raise if any DAG edge is violated (child > parent on present values).

    Edges are taken from :data:`_DAG_EDGES`. Equal values are allowed
    (non-increasing, not strictly decreasing). Edges where either endpoint
    is NaN for a sample are skipped for that sample.
    """
    for sample, row in df.iterrows():
        for parent, child in _DAG_EDGES:
            pv, cv = row[parent], row[child]
            if pd.notna(pv) and pd.notna(cv) and cv > pv:
                raise AssertionError(
                    f"summary_stats DAG violated for sample {sample!r}: "
                    f"{parent}={pv} < {child}={cv}"
                )


def _assert_strict_identities(
    df: pd.DataFrame,
    group_total: pd.Series | None,
    group_assigned: pd.Series | None,
    matrix_raw_rowsum: pd.Series | None,
) -> None:
    """Sanity checks enforced when ``strict=True``.

    Each check fires only when the corresponding source is provided AND
    the add_tags column is present for the sample. The checks express:

    * ``counted_depth`` (from add_tags) == row count in ``group.tsv.gz``
    * ``counted_assigned_depth`` (from add_tags) == row count in
      ``group.tsv.gz`` after the ``Unassigned,`` filter
    * ``counted_assigned_depth`` (from add_tags) == row sum of the
      pre-dedup count matrix
    """

    def _check(col: str, expected: pd.Series, label: str) -> None:
        observed = df[col]
        for sample in observed.index.intersection(expected.index):
            obs = observed.loc[sample]
            exp = expected.loc[sample]
            if pd.isna(obs) or pd.isna(exp):
                continue
            if int(obs) != int(exp):
                raise AssertionError(
                    f"summary_stats strict check failed for sample "
                    f"{sample!r}: {col}={int(obs)} from add_tags log "
                    f"!= {int(exp)} from {label}"
                )

    if group_total is not None:
        _check("counted_depth", group_total, "group.tsv.gz row count")
    if group_assigned is not None:
        _check(
            "counted_assigned_depth",
            group_assigned,
            "group.tsv.gz rows with non-Unassigned gene",
        )
    if matrix_raw_rowsum is not None:
        _check(
            "counted_assigned_depth",
            matrix_raw_rowsum,
            "row sum of count_matrix_raw",
        )


def summary_stats(
    *,
    fastq_stats: str | None = None,
    fastp_stats: str | None = None,
    star_logs: str | list[str] | None = None,
    add_tags_logs: str | list[str] | None = None,
    group_tsvs: str | list[str] | None = None,
    count_matrix_raw: str | None = None,
    count_matrix_dedup: str | None = None,
    strict: bool = False,
) -> pd.DataFrame:
    """Aggregate per-sample pipeline-stage counts into one DataFrame.

    Args:
        fastq_stats: seqkit ``--tabular`` stats file over raw input FASTQs
            (``foli qc`` input). Drives the ``raw_depth`` column.
        fastp_stats: seqkit ``--tabular`` stats file over fastp-trimmed FASTQs
            (``foli qc`` output). Drives the ``pass_qc_depth`` column.
        star_logs: Glob/path/list of STAR ``Log.final.out`` files. Each parent
            directory name is taken as the sample. Drives ``long_read_depth``
            (``Number of input reads``).
        add_tags_logs: Glob/path/list of ``foli_add_tags`` log files written
            with ``--log``. Samples come from the ``cell_tag`` field in each
            SUMMARY line. Drives ``not_na_adapter_depth``, ``good_umi_depth``,
            ``mapped_depth``, ``assigned_depth``, ``counted_depth``, and
            ``counted_assigned_depth``.
        group_tsvs: Glob/path/list of ``umi_tools group --group-out`` TSVs
            (``<sample>.group.tsv.gz``). Used only when ``strict=True`` to
            sanity-check ``counted_depth`` (row count) and
            ``counted_assigned_depth`` (rows after the ``Unassigned,``
            filter) against the values reported by ``foli_add_tags``.
        count_matrix_raw: Pre-UMI-dedup count matrix from
            ``foli get-count-mtx --output-raw``. When ``add_tags_logs`` is
            absent its row sums populate ``counted_assigned_depth`` as a
            backward-compatible fallback. When ``strict=True`` and both
            sources are present, the row sum is asserted equal to
            ``counted_assigned_depth`` from the add_tags log.
        count_matrix_dedup: UMI-dedup count matrix from
            ``foli get-count-mtx --output``. Row sums drive ``n_umi``; the
            per-row count of nonzero genes drives ``n_genes``.
        strict: If True, in addition to the always-on DAG check, assert
            that the ``counted`` / ``counted_assigned`` values reported by
            ``foli_add_tags`` agree with the ``group.tsv.gz`` row counts
            and the count-matrix row sums on every sample where both
            sources are present. Default False.

    Returns:
        DataFrame indexed by sample with columns in :data:`METRIC_COLUMNS`
        order. Dtype is pandas ``Int64`` so unset metrics are preserved as
        ``pd.NA``.

    Raises:
        AssertionError: If any DAG edge has child > parent on present
            values, or — when ``strict=True`` — if any of the cross-source
            identities described above is violated.
    """
    series_by_col: dict[str, pd.Series] = {}

    if fastq_stats is not None:
        series_by_col["raw_depth"] = _read_seqkit_r1_counts(fastq_stats)
    if fastp_stats is not None:
        series_by_col["pass_qc_depth"] = _read_seqkit_r1_counts(fastp_stats)
    if star_logs is not None:
        series_by_col["long_read_depth"] = _read_star_input_reads(
            _as_path_list(star_logs, suffix="Log.final.out")
        )
    if add_tags_logs is not None:
        parsed = _read_add_tags_summary(_as_path_list(add_tags_logs, suffix=None))
        for key in _ADD_TAGS_FIELDS:
            if not parsed[key].empty:
                series_by_col[f"{key}_depth"] = parsed[key]

    matrix_raw_rowsum: pd.Series | None = None
    if count_matrix_raw is not None:
        df_raw = _read_count_matrix(count_matrix_raw)
        matrix_raw_rowsum = df_raw.sum(axis=1).astype("int64")
        # Backward-compatible fallback: if add_tags didn't drive
        # counted_assigned_depth, populate it from the matrix row sum.
        if "counted_assigned_depth" not in series_by_col:
            series_by_col["counted_assigned_depth"] = matrix_raw_rowsum

    if count_matrix_dedup is not None:
        df_dedup = _read_count_matrix(count_matrix_dedup)
        series_by_col["n_umi"] = df_dedup.sum(axis=1).astype("int64")
        series_by_col["n_genes"] = (df_dedup > 0).sum(axis=1).astype("int64")

    samples = sorted({s for series in series_by_col.values() for s in series.index})
    df = pd.DataFrame(index=pd.Index(samples, name="sample"), columns=METRIC_COLUMNS)
    for col, series in series_by_col.items():
        df[col] = series.reindex(df.index)
    df = df.astype("Int64")

    _assert_dag_non_increasing(df)

    if strict:
        group_total = group_assigned = None
        if group_tsvs is not None:
            group_total, group_assigned = _read_group_tsv_counts(
                _as_path_list(group_tsvs, suffix=None)
            )
        _assert_strict_identities(
            df, group_total, group_assigned, matrix_raw_rowsum
        )

    return df
