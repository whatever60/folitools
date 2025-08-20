import subprocess
import re
import zipfile
from io import StringIO
from pathlib import Path
import tempfile
import os

import numpy as np
import pandas as pd
import polars as pl
from Bio import SeqIO


def extract_q30_from_fastqc_zip(zip_path: Path, inner_txt_path: str) -> float:
    """
    Extract the Q30 ratio from a FastQC zip archive.

    Args:
        zip_path: Path to the FastQC ZIP file.
        inner_txt_path: Path inside the ZIP to the `fastqc_data.txt` file,
            e.g. "{sample}_R1_001_fastqc/fastqc_data.txt".

    Returns:
        The fraction of bases with quality ≥ 30.
    """
    with zipfile.ZipFile(zip_path) as zf:
        with zf.open(inner_txt_path) as f:
            content = f.read().decode()

    match = re.search(
        r">>Per sequence quality scores.*?\n(.*?)>>END_MODULE",
        content,
        re.DOTALL,
    )
    if not match:
        raise ValueError(
            f"'Per sequence quality scores' section not found in {zip_path}"
        )

    section = match.group(1).strip()
    df = pd.read_csv(StringIO(section), sep="\t", comment=">")
    return df.loc[df["#Quality"] >= 30, "Count"].sum() / df["Count"].sum()


def read_fasta(fasta_path: Path) -> dict[str, str]:
    """
    Read a FASTA file into a dict of {record_id: sequence}.

    Args:
        fasta_path: Path to the FASTA file.

    Returns:
        Dictionary mapping sequence IDs to sequences.
    """
    return {r.id: str(r.seq) for r in SeqIO.parse(str(fasta_path), "fasta")}


def read_stat(path: Path | str) -> pd.DataFrame:
    """
    Read a simple two-column stats file and annotate sample and read.

    Args:
        path: Path to a tab-delimited file with an index column.

    Returns:
        DataFrame with added 'sample' and 'read' columns.
    """
    df = pd.read_table(path, index_col=0)
    df["sample"] = df.index.map(lambda x: Path(x).name.split(".")[0].split("_")[0])
    df["read"] = df.index.map(lambda x: "r1" if "R1_001" in x or "_1." in x else "r2")
    return df


def deduplicate_umi(df: pl.DataFrame) -> pl.DataFrame:
    """
    Count unique UMIs per gene.

    Args:
        df: Polars DataFrame with columns 'gene' and 'unique_id'.

    Returns:
        Polars DataFrame with columns ['gene', 'umi_count'].
    """
    return (
        df.select(["gene", "unique_id"])
        .unique()
        .group_by("gene")
        .len()
        .rename({"len": "umi_count"})
    )


def count_read_pairs(bam_path: str) -> int:
    """
    Count paired-end reads in a BAM/SAM by running samtools fastq once,
    saving R1/R2 to temporary FASTQ files, and counting them.

    Args:
        bam_path: Path to the input BAM or SAM file.

    Returns:
        The number of read pairs (asserting R1 count == R2 count).

    Raises:
        subprocess.CalledProcessError: If samtools fails.
        AssertionError: If the two counts disagree.
    """
    bam = Path(bam_path)
    if not bam.is_file():
        raise FileNotFoundError(f"BAM file not found: {bam_path}")

    # create two temp files for R1 and R2
    with (
        tempfile.NamedTemporaryFile(prefix="r1_", suffix=".fq", delete=False) as fwd,
        tempfile.NamedTemporaryFile(prefix="r2_", suffix=".fq", delete=False) as rev,
    ):
        r1_path = fwd.name
        r2_path = rev.name

    try:
        # run samtools fastq once, writing both mates
        subprocess.run(
            ["samtools", "fastq", "-1", r1_path, "-2", r2_path, bam_path],
            check=True,
            # stderr=subprocess.DEVNULL,
        )

        # count lines in each FASTQ and convert to reads
        def _count_reads(fastq_path: str) -> int:
            line_count = 0
            with open(fastq_path, "r") as fh:
                for _ in fh:
                    line_count += 1
            if line_count % 4 != 0:
                raise ValueError(f"FASTQ line count not multiple of 4: {fastq_path}")
            return line_count // 4

        r1_count = _count_reads(r1_path)
        r2_count = _count_reads(r2_path)

        assert r1_count == r2_count, f"Read pair mismatch: R1={r1_count}, R2={r2_count}"
        return r1_count

    finally:
        # clean up temp files
        for tmp in (r1_path, r2_path):
            try:
                os.remove(tmp)
            except OSError:
                pass


class FoliQC:
    """
    Compute QC metrics for a Foli sequencing sample.

    Args:
        data_dir: Root directory containing all stats and outputs.
        sample_name: Identifier for this sample.
        probe_set: Probe set name (used for primer lookup).
        primer_fwd: Path to forward-primer definitions (FASTA/CSV/TSV) or list.
        primer_rev: Path to reverse-primer definitions (FASTA/CSV/TSV) or list.
    """

    def __init__(
        self,
        data_dir: str | Path,
        sample_name: str,
        probe_set: str,
        primer_fwd: str | list[str],
        primer_rev: str | list[str],
    ) -> None:
        self.data_dir = Path(data_dir)
        self.sample_name = sample_name
        self.probe_set = probe_set

        # parse primers
        self.primer_fwds = self._load_primer_list(primer_fwd)
        self.primer_revs = self._load_primer_list(primer_rev)

        # initialize all metrics to None
        self.num_reads_raw: int | None = None
        self.num_reads_qc: int | None = None
        self.num_reads_primer: int | None = None
        self.num_reads_good_umi: int | None = None
        self.num_reads_mapped: int | None = None
        self.num_counts: int | None = None
        self.num_umi: int | None = None
        self.num_genes: int | None = None
        self.q30_r1: float | None = None
        self.q30_r2: float | None = None

        # load and compute
        self._load_stats()
        self._load_parquet_info()
        self._load_counts()
        self._compute_q30()
        self._compute_metrics()

    def _load_primer_list(self, src: str | list[str]) -> list[str]:
        """
        Parse a primer definition source into a list of primer IDs.

        Args:
            src: Path to FASTA/CSV/TSV file or a list of primer names.

        Returns:
            List of primer identifiers.
        """
        if isinstance(src, list):
            return src
        path = Path(src)
        suffix = path.suffix.lower().lstrip(".")
        if suffix in ("fasta", "fa", "fna"):
            return list(read_fasta(path).keys())
        elif suffix == "csv":
            return pd.read_csv(path).iloc[:, 0].astype(str).tolist()
        elif suffix in ("tsv", "txt"):
            return pd.read_table(path).iloc[:, 0].astype(str).tolist()
        else:
            raise ValueError(f"Unsupported primer file format: {suffix}")

    def _load_stats(self) -> None:
        """Read all *_stats files into attributes (or set to None)."""
        stat_files = {
            "fastq": "fastq.stats",
            "fastp": "fastp.stats",
            "rest": "rest.stats",
            "rest_all": "rest_all.stats",
        }
        for attr, fname in stat_files.items():
            path = self.data_dir / fname
            df = read_stat(path) if path.is_file() else None
            setattr(self, f"df_stat_{attr}", df)

        star_path = self.data_dir / "star" / self.sample_name / "ReadsPerGene.out.tab"
        if star_path.is_file():
            self.df_stat_star = pd.read_table(
                star_path,
                index_col=0,
                header=None,
                names=["unstranded", "forward", "reverse"],
            )
        else:
            self.df_stat_star = None

    def _load_parquet_info(self) -> None:
        """Read the rest_all_stats parquet and filter primer/UMI info."""
        p = self.data_dir / "rest_all_stats" / f"{self.sample_name}.parquet"
        if p.is_file():
            df = pl.read_parquet(p)
            self.df_info_rest_all = df

            # reads with any primer
            mask = (pl.col("primer_fwd") != "no_adapter") & (
                pl.col("primer_rev") != "no_adapter"
            )
            self.num_reads_primer = int(df.filter(mask).shape[0])

            # reads with good UMIs
            mask_good = (
                mask
                & (pl.col("umi5") != "")
                & (pl.col("umi3") != "")
                & (pl.col("umi5").str.len_chars() == 6)
                & (pl.col("umi3").str.len_chars() == 6)
                & (~pl.col("umi5").str.contains("N"))
                & (~pl.col("umi3").str.contains("N"))
            )
            df_good = df.filter(mask_good)
            self.df_info_rest_all_good = df_good
            self.num_reads_good_umi = int(df_good.shape[0])
        else:
            self.df_info_rest_all = None
            self.df_info_rest_all_good = None

    def _load_counts(self) -> None:
        """Read counts and deduplicate UMIs."""
        p = self.data_dir / "counts" / f"{self.sample_name}.group.tsv.gz"
        if p.is_file():
            df_counts = pl.read_csv(p, separator="\t")
            self.df_counts = df_counts
            self.num_counts = int(df_counts.shape[0])
            df_umi = deduplicate_umi(df_counts)
            self.num_umi = int(df_umi["umi_count"].sum())
            self.num_genes = int(df_umi.shape[0])
        else:
            self.df_counts = None

    def _compute_q30(self) -> None:
        """Compute Q30 for R1 and R2 from FastQC outputs."""
        fq_dir = self.data_dir / "fastq_fastqc"
        if fq_dir.is_dir():
            # R1
            r1_zips = list(fq_dir.glob(f"{self.sample_name}_*_R1_001_fastqc.zip"))
            if r1_zips:
                inner = f"{r1_zips[0].stem}/fastqc_data.txt"
                self.q30_r1 = extract_q30_from_fastqc_zip(r1_zips[0], inner)
            # R2
            r2_zips = list(fq_dir.glob(f"{self.sample_name}_*_R2_001_fastqc.zip"))
            if r2_zips:
                inner = f"{r2_zips[0].stem}/fastqc_data.txt"
                self.q30_r2 = extract_q30_from_fastqc_zip(r2_zips[0], inner)

    def _compute_metrics(self) -> None:
        """Compute read counts at each processing step."""
        # raw reads
        # df = self.df_stat_fastq
        df = getattr(self, "df_stat_fastq", None)
        if df is not None:
            cond = (df["sample"] == self.sample_name) & (df["read"] == "r1")
            if cond.any():
                self.num_reads_raw = int(df.loc[cond, "num_seqs"].item())

        # QC-passed reads
        # df = self.df_stat_fastp
        df = getattr(self, "df_stat_fastp", None)
        if df is not None:
            cond = (df["sample"] == self.sample_name) & (df["read"] == "r1")
            if cond.any():
                self.num_reads_qc = int(df.loc[cond, "num_seqs"].item())

        # mapped reads
        df = self.df_stat_star
        if df is not None:
            df2 = df.drop(index="N_unmapped", errors="ignore")
            self.num_reads_mapped = int(df2["unstranded"].sum())
