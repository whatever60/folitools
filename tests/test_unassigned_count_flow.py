"""Tests for unassigned BAM tagging and count-matrix filtering."""

from pathlib import Path
import tomllib

import pandas as pd
import pysam

from folitools import __version__
from folitools.add_tags import add_tags_wo_fastq
from folitools.cli import get_count_mtx
from folitools.get_matrix import read_counts


def test_add_tags_keeps_unassigned_primary_pair(tmp_path: Path) -> None:
    """Primary pairs with no assigned genes should still have R1 stamped with XF.

    umi_tools group/count in --paired mode only inspects R1, so add_tags
    writes CB/UC/XF to R1 only; R2 is passed through untouched.
    """
    input_bam = tmp_path / "input.bam"
    output_bam = tmp_path / "output.bam"

    header = {
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "SQ": [{"SN": "chr1", "LN": 1000}],
    }

    with pysam.AlignmentFile(input_bam, "wb", header=header) as bam_out:
        read1 = pysam.AlignedSegment()
        read1.query_name = "read1_AAAAAA_TTTTTT_FGR+FGR"
        read1.query_sequence = "A" * 50
        read1.flag = 99
        read1.reference_id = 0
        read1.reference_start = 100
        read1.mapping_quality = 60
        read1.cigarstring = "50M"
        read1.next_reference_id = 0
        read1.next_reference_start = 200
        read1.template_length = 150
        read1.query_qualities = pysam.qualitystring_to_array("I" * 50)
        read1.set_tag("XT", "Unassigned", value_type="Z")
        bam_out.write(read1)

        read2 = pysam.AlignedSegment()
        read2.query_name = "read1_AAAAAA_TTTTTT_FGR+FGR"
        read2.query_sequence = "T" * 50
        read2.flag = 147
        read2.reference_id = 0
        read2.reference_start = 200
        read2.mapping_quality = 60
        read2.cigarstring = "50M"
        read2.next_reference_id = 0
        read2.next_reference_start = 100
        read2.template_length = -150
        read2.query_qualities = pysam.qualitystring_to_array("I" * 50)
        read2.set_tag("XT", "Unassigned", value_type="Z")
        bam_out.write(read2)

    add_tags_wo_fastq(str(input_bam), str(output_bam), cell_tag="sample1")

    with pysam.AlignmentFile(output_bam, "rb") as bam_in:
        reads = list(bam_in.fetch(until_eof=True))

    assert len(reads) == 2
    r1 = next(r for r in reads if r.is_read1)
    r2 = next(r for r in reads if r.is_read2)

    assert r1.get_tag("XF") == "Unassigned,FGR+FGR"
    assert r1.get_tag("CB") == "sample1"
    assert r1.get_tag("UC") == "AAAAAATTTTTT"

    assert not r2.has_tag("XF")
    assert not r2.has_tag("CB")
    assert not r2.has_tag("UC")


def test_read_counts_filters_unassigned_gene_prefix(tmp_path: Path) -> None:
    """Rows whose gene starts with Unassigned should stay out of the matrix."""
    group_tsv = tmp_path / "sample.group.tsv"
    pd.DataFrame(
        {
            "read_id": ["read1", "read2"],
            "contig": ["chr1", "chr1"],
            "position": [100, 200],
            "gene": ["Unassigned,FGR+FGR", "GENE1,FGR+FGR"],
            "umi": ["AAAAAATTTTTT", "CCCCCCGGGGGG"],
            "umi_count": [1, 1],
            "final_umi": ["AAAAAATTTTTT", "CCCCCCGGGGGG"],
            "final_umi_count": [1, 1],
            "unique_id": ["0", "1"],
        }
    ).to_csv(group_tsv, sep="\t", index=False)

    matrix = read_counts([str(group_tsv)])

    assert list(matrix.columns) == ["GENE1"]
    assert matrix.iloc[0, 0] == 1


def test_read_counts_collapses_primer_annotations_before_counting(
    tmp_path: Path,
) -> None:
    """Primer variants for one gene should share raw and deduplicated counts."""
    group_tsv = tmp_path / "sample.group.tsv"
    pd.DataFrame(
        {
            "read_id": ["read1", "read2", "read3"],
            "contig": ["chr1", "chr1", "chr1"],
            "position": [100, 100, 100],
            "gene": [
                "GENE1,FGR1+RVR1",
                "GENE1,FGR2+RVR2",
                "GENE1,FGR1+RVR1",
            ],
            "umi": ["AAAA", "AAAA", "CCCC"],
            "umi_count": [2, 2, 1],
            "final_umi": ["AAAA", "AAAA", "CCCC"],
            "final_umi_count": [2, 2, 1],
            "unique_id": ["0", "0", "1"],
        }
    ).to_csv(group_tsv, sep="\t", index=False)

    raw = read_counts([str(group_tsv)], dedup_umi=False)
    dedup = read_counts([str(group_tsv)], dedup_umi=True)

    assert list(raw.columns) == ["GENE1"]
    assert raw.iloc[0, 0] == 3
    assert dedup.iloc[0, 0] == 2


def test_read_counts_sums_gene_ids_with_the_same_symbol(tmp_path: Path) -> None:
    """Distinct Ensembl IDs mapping to one symbol should preserve total reads."""
    group_tsv = tmp_path / "sample.group.tsv"
    gtf = tmp_path / "genes.gtf"
    pd.DataFrame(
        {
            "read_id": ["read1", "read2"],
            "contig": ["chr1", "chr1"],
            "position": [100, 200],
            "gene": [
                "ENSG00000000001.1,FGR1+RVR1",
                "ENSG00000000002.1,FGR2+RVR2",
            ],
            "umi": ["AAAA", "CCCC"],
            "umi_count": [1, 1],
            "final_umi": ["AAAA", "CCCC"],
            "final_umi_count": [1, 1],
            "unique_id": ["0", "1"],
        }
    ).to_csv(group_tsv, sep="\t", index=False)
    gtf.write_text(
        'chr1\ttest\tgene\t1\t100\t.\t+\t.\tgene_id "ENSG00000000001.1"; '
        'gene_name "GENE1";\n'
        'chr1\ttest\tgene\t200\t300\t.\t+\t.\tgene_id "ENSG00000000002.1"; '
        'gene_name "GENE1";\n'
    )

    matrix = read_counts([str(group_tsv)], gtf=str(gtf), dedup_umi=False)

    assert list(matrix.columns) == ["GENE1"]
    assert matrix.iloc[0, 0] == 2


def test_get_count_mtx_writes_package_version_header(tmp_path: Path) -> None:
    """The exported matrix should stamp the package version in the first cell."""
    group_tsv = tmp_path / "sample.group.tsv"
    output_tsv = tmp_path / "foli_counts.tsv"
    with (Path(__file__).resolve().parents[1] / "pyproject.toml").open("rb") as f:
        project_version = tomllib.load(f)["project"]["version"]
    pd.DataFrame(
        {
            "read_id": ["read1"],
            "contig": ["chr1"],
            "position": [100],
            "gene": ["GENE1,FGR+FGR"],
            "umi": ["AAAAAATTTTTT"],
            "umi_count": [1],
            "final_umi": ["AAAAAATTTTTT"],
            "final_umi_count": [1],
            "unique_id": ["0"],
        }
    ).to_csv(group_tsv, sep="\t", index=False)

    get_count_mtx(input_=[str(group_tsv)], output=str(output_tsv))

    assert __version__ == project_version
    assert output_tsv.read_text().splitlines()[0] == f"folitools {project_version}\tGENE1"
