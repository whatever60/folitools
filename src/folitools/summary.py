"""Per-sample pipeline-stage counts for foli runs.

``summary_stats`` aggregates read counts from each stage of the foli pipeline
into a single DataFrame (rows = samples, columns = metrics). Any source left
as ``None`` yields an all-NaN column. The columns are ordered so that, for a
healthy run, values are monotonically non-increasing across each row; an
``AssertionError`` is raised when that invariant is violated so pipeline
regressions surface immediately.
"""

from pathlib import Path
import re

import pandas as pd

from .utils import expand_path_to_list


METRIC_COLUMNS = (
    "raw_depth",
    "pass_qc_depth",
    "long_read_depth",
    "good_umi_depth",
    "not_na_adapter_depth",
    "properly_mapped_depth",
    "n_umi",
    "n_genes",
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
    """Parse STAR ``Log.final.out`` files for 'Number of input reads'.

    Sample name is inferred from the parent directory, matching the
    ``star_foli/<sample>/Log.final.out`` layout written by ``foli map``.
    """
    result: dict[str, int] = {}
    for p in paths:
        sample = Path(p).parent.name
        with open(p) as fh:
            for line in fh:
                m = _STAR_INPUT_RE.search(line)
                if m:
                    result[sample] = int(m.group(1))
                    break
    return pd.Series(result, dtype="int64")


def _read_add_tags_summary(paths: list[str]) -> tuple[pd.Series, pd.Series]:
    """Parse SUMMARY lines from ``foli_add_tags --log`` outputs."""
    good_umi: dict[str, int] = {}
    not_na_adapter: dict[str, int] = {}
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
                good_umi[sample] = int(fields["good_umi"])
                not_na_adapter[sample] = int(fields["not_na_adapter"])
    return (
        pd.Series(good_umi, dtype="int64"),
        pd.Series(not_na_adapter, dtype="int64"),
    )


def _read_count_matrix(path: str) -> pd.DataFrame:
    """Load a ``foli get-count-mtx`` matrix (samples x genes) into a DataFrame."""
    sep = "," if path.endswith((".csv", ".csv.gz")) else "\t"
    df = pd.read_csv(path, sep=sep, index_col=0)
    df.index.name = "sample"
    return df


def _as_path_list(arg: str | list[str] | None, suffix: str | None) -> list[str]:
    if arg is None:
        return []
    return expand_path_to_list(arg, suffix=suffix or "")


def _assert_monotonic_non_increasing(df: pd.DataFrame) -> None:
    """Raise if any row has a strictly-increasing step between present values.

    NaN entries are skipped; the comparison only looks at adjacent present
    values. Equal values are allowed (non-increasing, not strictly decreasing).
    """
    for sample, row in df.iterrows():
        values = [(col, v) for col, v in row.items() if pd.notna(v)]
        for (_, prev_v), (curr_c, curr_v) in zip(values, values[1:]):
            if curr_v > prev_v:
                raise AssertionError(
                    f"summary_stats not monotonically non-increasing for "
                    f"sample {sample!r}: {values} (violation at {curr_c})"
                )


def summary_stats(
    *,
    fastq_stats: str | None = None,
    fastp_stats: str | None = None,
    star_logs: str | list[str] | None = None,
    add_tags_logs: str | list[str] | None = None,
    count_matrix_raw: str | None = None,
    count_matrix_dedup: str | None = None,
) -> pd.DataFrame:
    """Aggregate per-sample pipeline-stage counts into one DataFrame.

    Args:
        fastq_stats: seqkit ``--tabular`` stats file over raw input FASTQs
            (``foli qc`` input). Drives the ``raw_depth`` column.
        fastp_stats: seqkit ``--tabular`` stats file over fastp-trimmed FASTQs
            (``foli qc`` output). Drives the ``pass_qc_depth`` column.
        star_logs: Glob/path/list of STAR ``Log.final.out`` files. Each parent
            directory name is taken as the sample. Drives ``long_read_depth``.
        add_tags_logs: Glob/path/list of ``foli_add_tags`` log files written
            with ``--log``. Samples come from the ``cell_tag`` field in each
            SUMMARY line. Drives ``good_umi_depth`` and ``not_na_adapter_depth``.
        count_matrix_raw: Pre-UMI-dedup count matrix from
            ``foli get-count-mtx --output-raw``. Row sums drive
            ``properly_mapped_depth``.
        count_matrix_dedup: UMI-dedup count matrix from
            ``foli get-count-mtx --output``. Row sums drive ``n_umi``; the
            per-row count of nonzero genes drives ``n_genes``.

    Returns:
        DataFrame indexed by sample with columns in :data:`METRIC_COLUMNS`
        order. Dtype is pandas ``Int64`` so unset metrics are preserved as
        ``pd.NA``.

    Raises:
        AssertionError: If any row violates monotonic non-increasing order
            across the present (non-NaN) metrics.
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
        good_umi, not_na_adapter = _read_add_tags_summary(
            _as_path_list(add_tags_logs, suffix=None)
        )
        series_by_col["good_umi_depth"] = good_umi
        series_by_col["not_na_adapter_depth"] = not_na_adapter
    if count_matrix_raw is not None:
        df_raw = _read_count_matrix(count_matrix_raw)
        series_by_col["properly_mapped_depth"] = df_raw.sum(axis=1).astype("int64")
    if count_matrix_dedup is not None:
        df_dedup = _read_count_matrix(count_matrix_dedup)
        series_by_col["n_umi"] = df_dedup.sum(axis=1).astype("int64")
        series_by_col["n_genes"] = (df_dedup > 0).sum(axis=1).astype("int64")

    samples = sorted({s for series in series_by_col.values() for s in series.index})
    df = pd.DataFrame(index=pd.Index(samples, name="sample"), columns=METRIC_COLUMNS)
    for col, series in series_by_col.items():
        df[col] = series.reindex(df.index)
    df = df.astype("Int64")

    _assert_monotonic_non_increasing(df)
    return df
