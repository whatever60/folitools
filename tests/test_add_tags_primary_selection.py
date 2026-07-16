"""Regression tests for non-compliant STAR primary-mate records."""

from io import StringIO
from pathlib import Path

import pysam
import pytest

from folitools._rust_native import run_add_tags
from folitools.add_tags import add_tags_wo_fastq


def _read(
    query_name: str,
    flag: int,
    position: int,
    *,
    xt: str = "Unassigned",
    hi: int | None = None,
) -> pysam.AlignedSegment:
    """Build one small mapped or unmapped alignment record."""
    read = pysam.AlignedSegment()
    read.query_name = query_name
    read.query_sequence = "A" * 50
    read.flag = flag
    read.reference_id = -1 if flag & 0x4 else 0
    read.reference_start = -1 if flag & 0x4 else position
    read.mapping_quality = 0 if flag & 0x4 else 60
    read.cigarstring = None if flag & 0x4 else "50M"
    read.next_reference_id = 0
    read.next_reference_start = 100
    read.template_length = 0
    read.query_qualities = pysam.qualitystring_to_array("I" * 50)
    read.set_tag("XT", xt, value_type="Z")
    if hi is not None:
        read.set_tag("HI", hi, value_type="i")
    return read


def _write_bam(path: Path, reads: list[pysam.AlignedSegment]) -> None:
    """Write records in the supplied order to an unsorted BAM."""
    header = {
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "SQ": [{"SN": "chr1", "LN": 1000}],
    }
    with pysam.AlignmentFile(path, "wb", header=header) as bam_out:
        for read in reads:
            bam_out.write(read)


def _run_successful_tagger(
    runner: str,
    input_bam: Path,
    output_bam: Path,
    log_path: Path,
) -> dict[str, int] | None:
    """Run one tagger and return Rust summary counters when available."""
    if runner == "rust":
        status = run_add_tags(
            [
                "foli_add_tags",
                "--input",
                str(input_bam),
                "--output",
                str(output_bam),
                "--cell_tag",
                "sample1",
                "--log",
                str(log_path),
            ]
        )
        assert status == 0
        summary = next(
            line for line in log_path.read_text().splitlines() if line.startswith("SUMMARY ")
        )
        return {
            key: int(value)
            for key, value in (field.split("=", 1) for field in summary.split()[2:])
        }

    log = StringIO()
    add_tags_wo_fastq(
        str(input_bam),
        str(output_bam),
        cell_tag="sample1",
        log_fh=log,
    )
    log_path.write_text(log.getvalue())
    return None


def _output_reads(path: Path) -> list[pysam.AlignedSegment]:
    """Read every output alignment without requiring an index."""
    with pysam.AlignmentFile(path, "rb") as bam_in:
        return list(bam_in.fetch(until_eof=True))


@pytest.mark.parametrize("runner", ["rust", "python"])
@pytest.mark.parametrize("ghost_first", [True, False])
def test_duplicate_r2_keeps_unique_mapped_candidate_regardless_of_order(
    tmp_path: Path,
    runner: str,
    ghost_first: bool,
) -> None:
    """STAR's unmapped R2 ghost must never displace the mapped primary R2."""
    query_name = "ghost-order_AAAAAA_TTTTTT_FGR+FGR"
    r1 = _read(query_name, 0x53, 100)
    ghost = _read(query_name, 0x85, -1, hi=0)
    mapped_r2 = _read(query_name, 0xA3, 200, hi=1)
    r2_records = [ghost, mapped_r2] if ghost_first else [mapped_r2, ghost]
    input_bam = tmp_path / f"{runner}.input.bam"
    output_bam = tmp_path / f"{runner}.output.bam"
    log_path = tmp_path / f"{runner}.log"
    _write_bam(input_bam, [r1, *r2_records])

    counters = _run_successful_tagger(runner, input_bam, output_bam, log_path)
    reads = _output_reads(output_bam)
    primary_r2 = [
        read
        for read in reads
        if read.is_read2 and not read.is_secondary and not read.is_supplementary
    ]
    secondary_r2 = [read for read in reads if read.is_read2 and read.is_secondary]

    assert len(primary_r2) == 1
    assert not primary_r2[0].is_unmapped
    assert primary_r2[0].reference_start == 200
    assert len(secondary_r2) == 1
    assert secondary_r2[0].is_unmapped
    assert not secondary_r2[0].has_tag("HI")
    assert "WARNING non-compliant primary count" in log_path.read_text()
    if counters is not None:
        assert counters["mapped"] == 1
        assert counters["counted"] == 1


@pytest.mark.parametrize("runner", ["rust", "python"])
def test_assigned_and_unassigned_pairs_have_expected_tags_and_counters(
    tmp_path: Path,
    runner: str,
) -> None:
    """Gene assignment should affect only assigned and counted-assigned totals."""
    assigned = "assigned_AAAAAA_TTTTTT_FGR+FGR"
    unassigned = "unassigned_CCCCCC_GGGGGG_FGR+FGR"
    input_bam = tmp_path / f"{runner}.input.bam"
    output_bam = tmp_path / f"{runner}.output.bam"
    log_path = tmp_path / f"{runner}.log"
    _write_bam(
        input_bam,
        [
            _read(assigned, 99, 100, xt="ENSG000001"),
            _read(assigned, 147, 200, xt="ENSG000001"),
            _read(unassigned, 99, 300),
            _read(unassigned, 147, 400),
        ],
    )

    counters = _run_successful_tagger(runner, input_bam, output_bam, log_path)
    r1_by_name = {
        read.query_name: read
        for read in _output_reads(output_bam)
        if read.is_read1 and not read.is_secondary and not read.is_supplementary
    }

    assert r1_by_name["assigned"].get_tag("XF") == "ENSG000001,FGR+FGR"
    assert r1_by_name["unassigned"].get_tag("XF") == "Unassigned,FGR+FGR"
    if counters is not None:
        assert counters == {
            "total_r1": 2,
            "not_na_adapter": 2,
            "good_umi": 2,
            "mapped": 2,
            "assigned": 1,
            "counted": 2,
            "counted_assigned": 1,
        }


@pytest.mark.parametrize("runner", ["rust", "python"])
def test_duplicate_r1_keeps_mapped_candidate_and_counts_once(
    tmp_path: Path,
    runner: str,
) -> None:
    """A duplicate R1 keeps the mapped record and Rust counts the QNAME once."""
    query_name = "duplicate-r1_AAAAAA_TTTTTT_FGR+FGR"
    input_bam = tmp_path / f"{runner}.input.bam"
    output_bam = tmp_path / f"{runner}.output.bam"
    log_path = tmp_path / f"{runner}.log"
    _write_bam(
        input_bam,
        [
            _read(query_name, 0x45, -1, hi=0),
            _read(query_name, 99, 100, xt="ENSG000001", hi=1),
            _read(query_name, 147, 200, xt="ENSG000001", hi=1),
        ],
    )

    counters = _run_successful_tagger(runner, input_bam, output_bam, log_path)
    primary_r1 = [
        read
        for read in _output_reads(output_bam)
        if read.is_read1 and not read.is_secondary and not read.is_supplementary
    ]

    assert len(primary_r1) == 1
    assert not primary_r1[0].is_unmapped
    assert primary_r1[0].reference_start == 100
    if counters is not None:
        assert counters == {
            "total_r1": 1,
            "not_na_adapter": 1,
            "good_umi": 1,
            "mapped": 1,
            "assigned": 1,
            "counted": 1,
            "counted_assigned": 1,
        }


@pytest.mark.parametrize("runner", ["rust", "python"])
def test_ordinary_primary_pair_is_unchanged(tmp_path: Path, runner: str) -> None:
    """A compliant mapped pair should keep its flags and produce no warning."""
    query_name = "ordinary_AAAAAA_TTTTTT_FGR+FGR"
    input_bam = tmp_path / f"{runner}.input.bam"
    output_bam = tmp_path / f"{runner}.output.bam"
    log_path = tmp_path / f"{runner}.log"
    _write_bam(
        input_bam,
        [
            _read(query_name, 99, 100, xt="ENSG000001"),
            _read(query_name, 147, 200, xt="ENSG000001"),
        ],
    )

    counters = _run_successful_tagger(runner, input_bam, output_bam, log_path)
    reads = _output_reads(output_bam)

    assert sorted(read.flag for read in reads) == [99, 147]
    assert "WARNING" not in log_path.read_text()
    if counters is not None:
        assert counters["total_r1"] == 1
        assert counters["mapped"] == 1
        assert counters["counted"] == 1


@pytest.mark.parametrize("runner", ["rust", "python"])
@pytest.mark.parametrize(
    ("candidate_flags", "mapped_count"),
    [
        ((0xA3, 0xA3), 2),
        ((0x85, 0x85), 0),
    ],
)
def test_indistinguishable_primary_candidates_fail_as_ambiguous(
    tmp_path: Path,
    runner: str,
    candidate_flags: tuple[int, int],
    mapped_count: int,
) -> None:
    """Candidates without one unique mapped mate must fail instead of guessing."""
    query_name = "ambiguous_AAAAAA_TTTTTT_FGR+FGR"
    input_bam = tmp_path / f"{runner}.input.bam"
    output_bam = tmp_path / f"{runner}.output.bam"
    log_path = tmp_path / f"{runner}.log"
    positions = [200, 300]
    _write_bam(
        input_bam,
        [
            _read(query_name, 0x53, 100),
            *[
                _read(query_name, flag, -1 if flag & 0x4 else position)
                for flag, position in zip(candidate_flags, positions, strict=True)
            ],
        ],
    )

    if runner == "rust":
        status = run_add_tags(
            [
                "foli_add_tags",
                "--input",
                str(input_bam),
                "--output",
                str(output_bam),
                "--log",
                str(log_path),
            ]
        )
        assert status == 1
    else:
        log = StringIO()
        with pytest.raises(ValueError) as error:
            add_tags_wo_fastq(
                str(input_bam),
                str(output_bam),
                log_fh=log,
            )
        log_path.write_text(log.getvalue())
        message = str(error.value)
        assert "ambiguous" in message.lower()
        assert "R2" in message
        assert f"{mapped_count} mapped" in message
        assert query_name in message

    assert "WARNING non-compliant primary count" in log_path.read_text()
