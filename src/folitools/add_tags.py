import sys
from typing import Annotated

import pysam
from cyclopts import run, Parameter


def add_tags_wo_fastq(
    bam_input: str,
    bam_output: str,
    umi_tag_s_name: str = "US",  # single-read UMI
    umi_tag_s_correct_name: str = "UC",  # single-read UMI (correct ones only)
    cell_tag_name: str = "CB",
    primer_tag_name: str = "PR",
    gene_with_primer_tag_name: str = "XF",
    cell_tag: str | None = None,
    umi_length: int = 6,
) -> None:
    """
    Reads SAM from stdin, parses embedded UMI sequences from the read names,
    adds UMI and optional cell tags, and writes a BAM file.

    VERY IMPORTANT NOTE:
        The read IDs must be formatted as:
            <original_read_id>_{umi_r1}_{umi_r2}
        where the two UMI sequences are appended to the original read name,
        separated by underscores. These UMIs will be split out and tagged
        appropriately.

    Args:
        bam_input: Standard input stream for SAM data.
        bam_output: Path to the output BAM file ('-' or None for stdout).
        umi_tag_s_name: Name of the single-read UMI tag (default: 'US').
        umi_tag_p_name: Name of the paired-read UMI tag (default: 'UP').
        cell_tag_name: Name of the cell barcode tag (default: 'CB').
        adapter_tag_name: Name of the adapter tag (unused here) (default: 'XA').
        cell_tag: Optional cell barcode value to add to every read.
        umi_length: Expected length of each individual UMI (default: 6).
    """
    if bam_input == "-":
        bam_input = "/dev/stdin"
    mode_out = "w"
    if bam_output == "-":
        mode_out = "wb0"  # same as wbu, i.e., uncompressed BAM
    elif bam_output.endswith(".bam"):
        mode_out = "wb"
    elif bam_output.endswith(".sam"):
        mode_out = "w"
    else:
        raise ValueError("Output file must end with .bam or .sam")

    def process_and_write_pair():
        if len(primary_alignments) != 2:
            msg = f"Expected 2 primary alignments, got {len(primary_alignments)}"
            msg += f" (current read: {read_id}, next read: {query_name_current})"
            raise ValueError(msg)

        # Join unique XS tags, filter out any "Unassigned" tags if others exist
        genes = {
            gene for tag in xt_tags if tag != "Unassigned" for gene in tag.split(",")
        }
        if genes:
            gene_with_primer = ",".join(sorted(genes) + [primers])
            for r, read in primary_alignments:
                read.set_tag(
                    gene_with_primer_tag_name, gene_with_primer, value_type="Z"
                )
                bam_out.write(read)

    with (
        # NOTE: pysam read can handle both BAM and SAM format automatically if opened with "rb".
        pysam.AlignmentFile(bam_input, "rb") as sam_in,
        pysam.AlignmentFile(bam_output, mode_out, header=sam_in.header) as bam_out,
    ):
        primary_alignments: list[tuple[str, pysam.AlignedSegment]] = []
        xt_tags = []
        read_id = None
        for row_idx, read in enumerate(sam_in.fetch(until_eof=True)):
            query_name_current = read.query_name
            assert query_name_current is not None, "Missing query name in read"
            if read_id is not None and read_id != query_name_current:
                # We've encountered a new read ID, so we need to process the previous one.
                process_and_write_pair()
                primary_alignments = []
                xt_tags = []
                read_id = None

            if read_id is None:
                read_id = query_name_current

            # Add the cell tag
            if cell_tag:
                read.set_tag(cell_tag_name, cell_tag, value_type="Z")

            # if "_" not in query_name:
            #     # UMI sequences are not embedded in read name. We assume that these
            #     # reads already have their UMI extracted and put in the tag.
            #     if not read.has_tag(umi_tag_s_name):
            #         raise ValueError(f"Missing UMI tag in read {query_name}")
            #     bam_out.write(read)
            #     continue

            # Add the UMI tag from the FASTQ sequence
            query_name_simple, umi1, umi2, primers = query_name_current.split("_", 3)
            primer_fwd, primer_rev = primers.split("+")
            read.query_name = query_name_simple
            criteria = [
                len(umi1) == umi_length,
                len(umi2) == umi_length,
                "N" not in umi1,
                "N" not in umi2,
                primer_fwd != "no_adapter",
                primer_rev != "no_adapter",
                # primer_fwd == primer_rev,
            ]
            # if not all(criteria):
            #     seq1_c = seq2_c = "AAAAAA"
            # else:
            #     seq1_c, seq2_c = seq1, seq2
            if read.is_read1:
                read.set_tag(umi_tag_s_name, umi1, value_type="Z")
            else:
                read.set_tag(umi_tag_s_name, umi2, value_type="Z")

            read.set_tag(primer_tag_name, primers, value_type="Z")

            if all(criteria):
                read.set_tag(umi_tag_s_correct_name, umi1 + umi2, value_type="Z")

            if not read.has_tag("XT"):  # No feature assigned by featurecounts
                assert str(read.get_tag("XS")).startswith("Unassigned")
                read.set_tag("XN", -1, value_type="i")
                read.set_tag("XT", "Unassigned", value_type="Z")
            if read.is_secondary or read.is_supplementary:
                # Immediate write for non-primary alignments
                bam_out.write(read)
            else:
                key = "R1" if read.is_read1 else "R2"
                primary_alignments.append((key, read))
            xt_tags.append(str(read.get_tag("XT")))

        # Wrap up the last pair
        process_and_write_pair()


def main(
    input_: Annotated[
        str,
        Parameter(
            name=["--input", "-i"],
            help="Path to the input BAM file ('-' or omit for stdin)",
        ),
    ] = "-",
    output: Annotated[
        str,
        Parameter(
            name=["--output", "-o"],
            help="Path to the output BAM file ('-' or omit for stdout)",
        ),
    ] = "-",
    cell_tag_name: Annotated[
        str,
        Parameter(
            name="--cell_tag_name",
            help="Name of the cell tag to add (default: 'CB')",
        ),
    ] = "CB",
    cell_tag: Annotated[
        str | None,
        Parameter(
            name="--cell_tag",
            help="Value of the cell tag to add (default: None)",
        ),
    ] = None,
) -> int:
    """
    Dispatch to add_tags or add_tags_wo_fastq based on presence of --index_fastq.

    All existing comments and logic are preserved.
    """
    try:
        add_tags_wo_fastq(
            bam_input=input_,
            bam_output=output,
            cell_tag_name=cell_tag_name,
            cell_tag=cell_tag,
        )
    except Exception as e:
        print(e, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(run(main))
