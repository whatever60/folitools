"""Tests for sequencing filename and read-ID parsing helpers."""

import pytest

from folitools.read_names import (
    read_number_from_path,
    sample_name_from_path,
    split_tagged_qname,
    split_umi_read_id,
    strip_mate_suffix,
)


@pytest.mark.parametrize(
    "path,expected",
    [
        ("sample_S1_L001_R1_001.fastq.gz", "sample"),
        ("sample_S1_R2_001.fastq.gz", "sample"),
        ("sample_R1.fastq.gz", "sample"),
        ("sample_1.fq.gz", "sample"),
        ("patient_A.1.fq", "patient_A"),
        ("treatment_group_001_R1.fastq.gz", "treatment_group_001"),
        ("yq-foli-na-ctrl-batch2-YQ1_R1.fastq.gz", "yq-foli-na-ctrl-batch2-YQ1"),
        ("result.sorted.bam", "result"),
        ("complex_sample_name.sam", "complex_sample_name"),
    ],
)
def test_sample_name_from_path_removes_only_trailing_sequencing_tokens(
    path: str, expected: str
) -> None:
    """Sample extraction should keep meaningful underscores in sample IDs."""
    assert sample_name_from_path(path) == expected


@pytest.mark.parametrize(
    "path,expected",
    [
        ("sample_S1_L001_R1_001.fastq.gz", "r1"),
        ("sample_S1_L001_R2_001.fastq.gz", "r2"),
        ("sample_R1.fastq.gz", "r1"),
        ("sample_R2.fastq.gz", "r2"),
        ("sample_1.fq.gz", "r1"),
        ("sample_2.fq.gz", "r2"),
    ],
)
def test_read_number_from_path_accepts_illumina_and_aviti_names(
    path: str, expected: str
) -> None:
    """Read-number parsing should accept both Illumina and AVITI exports."""
    assert read_number_from_path(path) == expected


def test_split_umi_read_id_allows_underscores_in_original_id() -> None:
    """Post-cutadapt read IDs should split UMIs from the right."""
    assert split_umi_read_id("READ_with_underscores_ACGTAA_TTGGCC") == (
        "READ_with_underscores",
        "ACGTAA",
        "TTGGCC",
    )


def test_split_tagged_qname_allows_underscores_in_original_id() -> None:
    """Mapping QNAME parsing should split UMI and primer fields from the right."""
    assert split_tagged_qname("READ_with_underscores_ACGTAA_TTGGCC_FWD+REV") == (
        "READ_with_underscores",
        "ACGTAA",
        "TTGGCC",
        "FWD+REV",
    )


def test_strip_mate_suffix_removes_legacy_read_suffix() -> None:
    """Legacy `/1` and `/2` suffixes should not make mates look unrelated."""
    assert strip_mate_suffix("@READ/1") == "@READ"
    assert strip_mate_suffix("@READ/2") == "@READ"
    assert strip_mate_suffix("@READ") == "@READ"
