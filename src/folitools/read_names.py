"""Helpers for FASTQ/BAM sample names and Foli-annotated read IDs."""

from pathlib import Path
import base64
import re
import sys

import dnaio


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
_DNA_DIGITS = str.maketrans("ACGTN_", "012345")


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
    if "_~" in qname:
        read_id, suffix = qname.rsplit("_~", 1)
        encoded, primers = suffix.split("_", 1)
        packed = int.from_bytes(base64.urlsafe_b64decode(encoded.replace(".", "_") + "=" * (-len(encoded) % 4)), "big")
        bases = []
        while packed > 1:
            bases.append("ACGTN_"[packed & 7])
            packed >>= 3
        if packed != 1:
            raise ValueError(f"Invalid packed QNAME: {qname}")
        umi1, umi2 = "".join(reversed(bases)).split("_")
        return read_id, umi1, umi2, primers
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


def bam_safe_qname(qname: str) -> str:
    """Pack only overlong UMIs using three bits per A/C/G/T/N and separator."""
    if len(qname.encode()) <= 254:
        return qname
    read_id, umi1, umi2, primers = split_tagged_qname(qname)
    packed = int("1" + (umi1 + "_" + umi2).translate(_DNA_DIGITS), 8)
    encoded = base64.urlsafe_b64encode(packed.to_bytes((packed.bit_length() + 7) // 8, "big")).decode().rstrip("=").replace("_", ".")
    result = f"{read_id}_~{encoded}_{primers}"
    if len(result.encode()) > 254:
        raise ValueError("Read ID and primer names exceed BAM's QNAME capacity after UMI packing")
    return result


if __name__ == "__main__":
    with dnaio.open(sys.stdin.buffer, mode="r", interleaved=True) as source, dnaio.open(sys.argv[1], sys.argv[2], mode="w", fileformat="fastq") as output:
        for r1, r2 in source:
            r1.name = bam_safe_qname(r1.name.replace(" ", "_"))
            r2.name = bam_safe_qname(r2.name.replace(" ", "_"))
            output.write(r1, r2)
