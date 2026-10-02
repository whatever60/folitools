"""Lossless three-bit transport for overlong annotated read names."""

from pathlib import Path
import random
import subprocess
import sys

import dnaio
import pytest

from folitools.read_names import bam_safe_qname, split_tagged_qname


def test_three_bit_roundtrip_all_padding_lengths():
    """Preserve every base, leading A, Ns, and the two variable UMI lengths."""
    rng = random.Random(7)
    for length in range(151):
        umi1 = "".join(rng.choices("ACGTN", k=length))
        umi2 = "".join(rng.choices("ACGTN", k=150))
        expected = ("read_with_underscores", umi1, umi2, "FGR|other+primer_with_underscores")
        original = "_".join(expected)
        encoded = bam_safe_qname(original)
        assert len(encoded.encode()) <= 254
        assert split_tagged_qname(encoded) == expected
        assert "/" not in encoded
        if len(original) <= 254:
            assert encoded == original


def test_fastq_stream_changes_headers_only(tmp_path):
    """Keep paired sequences and qualities identical while packing long headers."""
    read_id = "AV100007:20260928_yq_mayo_batch2_ot_moa_training_ot_takeda_batch1:2613420416:2:11401:1359:3162"
    name = f"{read_id}_{'A' * 90}_{'N' * 90} XIAP+IL1R1"
    source = tmp_path / "input.fq"
    source.write_text(f"@{name}\nACGT\n+\nIJKL\n@{name}\nTGCA\n+\nLMNO\n")
    r1, r2 = tmp_path / "r1.fq", tmp_path / "r2.fq"
    with source.open("rb") as input_stream:
        subprocess.run([sys.executable, "-m", "folitools.read_names", str(r1), str(r2)], stdin=input_stream, check=True)
    with dnaio.open(r1, r2) as output:
        first, second = next(iter(output))
    assert (first.sequence, first.qualities) == ("ACGT", "IJKL")
    assert (second.sequence, second.qualities) == ("TGCA", "LMNO")
    assert first.name == second.name
    assert split_tagged_qname(first.name) == (read_id, "A" * 90, "N" * 90, "XIAP+IL1R1")


def test_unrepresentable_names_fail_before_writing_bam():
    """Never silently truncate names that cannot fit even with packed UMIs."""
    with pytest.raises(ValueError, match="QNAME capacity"):
        bam_safe_qname("_".join(["r" * 254, "A" * 90, "N" * 90, "F+R"]))
