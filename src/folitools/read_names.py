"""Helpers for FASTQ/BAM sample names and Foli-annotated read IDs."""

from pathlib import Path
import re


_PATH_SUFFIXES = (
    ".fastq.gz",
    ".fq.gz",
    ".fastq",
    ".fq",
    ".sam.gz",
    ".sam.bz2",
    ".bam",
    ".sam",
)
_READ_RE = re.compile(
    r"^(?P<sample>.+?)(?:_S\d+)?(?:_L\d{3})?[._](?P<read>R?[12])(?:_\d{3})?_*$",
    re.IGNORECASE,
)
QNAME_PRIMER_SEPARATOR = "|"


def strip_path_suffix(path: str | Path) -> str:
    """Return the basename after removing one supported sequencing-file suffix."""
    name = Path(path).name
    for suffix in _PATH_SUFFIXES:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def sample_name_from_path(path: str | Path) -> str:
    """Infer the sample name by removing only trailing read/lane/sort tokens."""
    name = strip_path_suffix(path)
    if name.endswith(".sorted"):
        name = name[: -len(".sorted")]
    match = _READ_RE.match(name)
    if match is not None:
        return match["sample"]
    return name


def read_number_from_path(path: str | Path) -> str:
    """Return ``r1`` or ``r2`` from a FASTQ filename."""
    name = strip_path_suffix(path)
    match = _READ_RE.match(name)
    if match is None:
        raise ValueError(f"Cannot infer read number from file name: {path}")
    read = match["read"].lower()
    if read.endswith("1"):
        return "r1"
    return "r2"


def strip_mate_suffix(read_id: str) -> str:
    """Remove legacy ``/1`` or ``/2`` mate suffixes from a FASTQ read ID."""
    if read_id.endswith(("/1", "/2")):
        return read_id[:-2]
    return read_id


def split_umi_read_id(read_id: str) -> tuple[str, str, str]:
    """Split ``<read_id>_<umi5>_<umi3>`` while allowing underscores in read_id."""
    parts = read_id.rsplit("_", 2)
    if len(parts) != 3:
        raise ValueError(f"Unexpected read ID format: {read_id}")
    return parts[0], parts[1], parts[2]


def split_tagged_qname(qname: str) -> tuple[str, str, str, str]:
    """Split ``<read_id>_<umi5>_<umi3>|<primer5+primer3>`` from a QNAME."""
    parts = qname.rsplit(QNAME_PRIMER_SEPARATOR, 1)
    if len(parts) != 2:
        raise ValueError(f"Unexpected read ID format: {qname}")
    read_id_with_umis, primers = parts
    if "+" not in primers:
        raise ValueError(f"Unexpected primer format in read ID: {qname}")
    read_id, umi1, umi2 = split_umi_read_id(read_id_with_umis)
    return read_id, umi1, umi2, primers
