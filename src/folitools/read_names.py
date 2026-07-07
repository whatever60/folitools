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
_UMI_RE = re.compile(r"^[ACGTN]*$", re.IGNORECASE)


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
    """Split ``<read_id>_<umi1>_<umi2>`` while allowing underscores in read_id."""
    parts = read_id.rsplit("_", 2)
    if len(parts) != 3:
        raise ValueError(f"Unexpected read ID format: {read_id}")
    return parts[0], parts[1], parts[2]


def split_tagged_qname(qname: str) -> tuple[str, str, str, str]:
    """Split the 0.7-style ``<read_id>_<umi1>_<umi2>_<primer1+primer2>`` QNAME."""
    parts = qname.split("_")
    for primer_idx in range(3, len(parts)):
        umi1 = parts[primer_idx - 2]
        umi2 = parts[primer_idx - 1]
        primers = "_".join(parts[primer_idx:])
        if (
            _UMI_RE.match(umi1) is not None
            and _UMI_RE.match(umi2) is not None
            and "+" in primers
        ):
            read_id = "_".join(parts[: primer_idx - 2])
            if read_id == "":
                raise ValueError(f"Unexpected read ID format: {qname}")
            return read_id, umi1, umi2, primers
    raise ValueError(f"Unexpected read ID format: {qname}")
