"""Tests for unassigned BAM tagging and count-matrix filtering."""

from pathlib import Path

import pandas as pd
import pysam

from folitools.add_tags import add_tags_wo_fastq
from folitools.get_matrix import read_counts


def test_add_tags_keeps_unassigned_primary_pair(tmp_path: Path) -> None:
    """Primary pairs with no assigned genes should still be written with XF."""
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
    assert all(read.get_tag("XF") == "Unassigned,FGR+FGR" for read in reads)
    assert all(read.get_tag("CB") == "sample1" for read in reads)
    assert all(read.get_tag("UC") == "AAAAAATTTTTT" for read in reads)


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
