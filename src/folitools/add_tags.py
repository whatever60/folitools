import sys
import gzip
from typing import Annotated, TextIO

import pysam
from Bio import SeqIO
from tqdm.auto import tqdm
from cyclopts import App, Parameter


app = App()


def _get_fastq_id_seq_interleaved(
    fastq_iter,
) -> tuple[str, str, str] | tuple[None, str, str]:
    fastq_record1 = next(fastq_iter, None)
    fastq_record2 = next(fastq_iter, None)
    if fastq_record1 is None:
        return None, "", ""
    else:
        if fastq_record2 is None:
            raise ValueError("Second FASTQ record is missing")
        if not fastq_record1.id == fastq_record2.id:
            raise ValueError(
                f"Interleaved FASTQ records do not match: {fastq_record1.id} != {fastq_record2.id}"
            )
        return fastq_record1.id, str(fastq_record1.seq), str(fastq_record2.seq)


def add_tags(
    sam_stream: TextIO,
    fastq_file: str,
    output_bam: str,
    umi_tag_s_name: str = "US",  # single-read UMI
    umi_tag_p_name: str = "UP",  # concatenated UMI from paired reads
    cell_tag_name: str = "CB",
    adapter_tag_name: str = "XA",
    cell_tag: str | None = None,
    umi_length: int | None = 6,
) -> None:
    """
    Reads SAM from stdin, adds UMI and optional cell tags from the FASTQ file, and writes a BAM file.

    VERY IMPORTANT NOTE:
        1. The sam stream comes from STAR and id of the reads has been modified to
            <original_read_id>+<adapter_name>.
        2. The FASTQ file is interleaved, with read id unmodified.
        3. Both the fastq inputs to STAR and the FASTQ here have been sorted by read id.
            And all read ids (before modification) in the STAR input are guaranteed to
            appear in the FASTQ file here.

    Args:
        sam_stream: Standard input stream for SAM data.
        fastq_file (str): Path to the second FASTQ file.
        output_bam (str): Path to the output BAM file with tags.
        umi_tag_name (str): Name of the UMI tag to add (default: 'UB').
        cell_tag_name (str): Name of the cell tag to add (default: 'CB').
        cell_tag (str): Value of the cell tag to add (default: 'default_cell_tag').
    """

    # Parse SAM from stdin and create a BAM file
    with (
        gzip.open(fastq_file, "rt") as f,
        pysam.AlignmentFile(sam_stream, "r") as sam_in,
        pysam.AlignmentFile(output_bam, "wb", header=sam_in.header) as bam_out,
    ):
        fastq_iter = iter(tqdm(SeqIO.parse(f, "fastq")))

        # Initialize a FASTQ record
        index_name, seq1, seq2 = _get_fastq_id_seq_interleaved(fastq_iter)
        for read in sam_in:
            read_name, adapter = read.query_name.split("+")
            read.query_name = read_name
            # read_name = read.query_name
            while read_name != index_name:
                index_name, seq1, seq2 = _get_fastq_id_seq_interleaved(fastq_iter)
                if index_name is None:
                    raise ValueError(
                        f"No more FASTQ records, current read: {read_name}"
                    )

            # read.set_tag(adapter_tag_name, adapter, value_type="Z")

            if read.flag & 0x40:  # read is first in pair
                read.set_tag(umi_tag_s_name, seq1, value_type="Z")
            else:
                read.set_tag(umi_tag_s_name, seq2, value_type="Z")
            # Add the UMI tag from the FASTQ sequence
            if umi_length is not None:
                criteria = [
                    len(seq1) == umi_length,
                    len(seq2) == umi_length,
                    "N" not in seq1,
                    "N" not in seq2,
                ]
                if all(criteria):
                    read.set_tag(umi_tag_p_name, seq1 + seq2, value_type="Z")

            # index_name_next, seq_next = _get_fastq_id_seq(fastq_iter)

            # for read in sam_in:
            #     read_name = read.query_name
            #     name_match = False
            #     if read_name == index_name:
            #         name_match = True
            #     elif read_name == index_name_next:
            #         name_match = True
            #         # advance index name
            #         index_name, seq = index_name_next, seq_next
            #         index_name_next, seq_next = _get_fastq_id_seq(fastq_iter)
            #     else:  # this read does not have UMI
            #         pass

            #     if name_match:
            #         # Add the UMI tag from the FASTQ sequence
            #         if umi_length is not None:
            #             if umi_length == len(seq):
            #                 read.set_tag(umi_tag_name, seq, value_type="Z")
            #             else:
            #                 read.set_tag(umi_tag_bad_name, seq, value_type="Z")

            # Add the cell tag
            if cell_tag:
                read.set_tag(cell_tag_name, cell_tag, value_type="Z")

            bam_out.write(read)
        # # Check if the FASTQ file still has unmatched entries
        # if index_name_next is not None:
        #     print(index_name, read.query_name)
        #     raise ValueError("Extra index sequences")


def add_tags_wo_fastq(
    sam_stream: TextIO,
    output_bam: str,
    umi_tag_s_name: str = "US",  # single-read UMI
    umi_tag_p_name: str = "UP",  # concatenated UMI from paired reads
    cell_tag_name: str = "CB",
    adapter_tag_name: str = "XA",
    cell_tag: str | None = None,
    umi_length: int | None = 6,
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
        sam_stream: Standard input stream for SAM data.
        output_bam: Path to the output BAM file ('-' or None for stdout).
        umi_tag_s_name: Name of the single-read UMI tag (default: 'US').
        umi_tag_p_name: Name of the paired-read UMI tag (default: 'UP').
        cell_tag_name: Name of the cell barcode tag (default: 'CB').
        adapter_tag_name: Name of the adapter tag (unused here) (default: 'XA').
        cell_tag: Optional cell barcode value to add to every read.
        umi_length: Expected length of each individual UMI (default: 6).
    """
    # open the SAM input
    sam_in = pysam.AlignmentFile(sam_stream, "r")
    # decide BAM output
    if output_bam in (None, "-"):
        bam_out = pysam.AlignmentFile("-", "wb0", header=sam_in.header)
    else:
        bam_out = pysam.AlignmentFile(output_bam, "wb", header=sam_in.header)

    with sam_in, bam_out:
        for read in sam_in:
            read_name, seq1, seq2 = read.query_name.split("_")
            if read.is_read1:
                read.set_tag(umi_tag_s_name, seq1, value_type="Z")
            else:
                read.set_tag(umi_tag_s_name, seq2, value_type="Z")
            read.query_name = read_name
            # Add the UMI tag from the FASTQ sequence
            if umi_length is not None:
                criteria = [
                    len(seq1) == umi_length,
                    len(seq2) == umi_length,
                    "N" not in seq1,
                    "N" not in seq2,
                ]
                if all(criteria):
                    read.set_tag(umi_tag_p_name, seq1 + seq2, value_type="Z")
            # Add the cell tag
            if cell_tag:
                read.set_tag(cell_tag_name, cell_tag, value_type="Z")
            bam_out.write(read)


@app.command(help="Add UMI and cell tags to BAM file from FASTQ sequences.")
def main(
    index_fastq: Annotated[
        str | None,
        Parameter(
            name="--index_fastq",
            help="Path to the interleaved FASTQ file (omit to use --add-tags-wo-fastq mode)",
        ),
    ] = None,
    output: Annotated[
        str | None,
        Parameter(
            name=["--output", "-o"],
            help="Path to the output BAM file ('-' or omit for stdout)",
        ),
    ] = None,
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
) -> None:
    """
    Dispatch to add_tags or add_tags_wo_fastq based on presence of --index_fastq.

    All existing comments and logic are preserved.
    """
    if index_fastq is None:
        add_tags_wo_fastq(
            sam_stream=sys.stdin,
            output_bam=output or "-",
            cell_tag_name=cell_tag_name,
            cell_tag=cell_tag,
        )
    else:
        add_tags(
            sam_stream=sys.stdin,
            fastq_file=index_fastq,
            output_bam=output or "-",
            cell_tag_name=cell_tag_name,
            cell_tag=cell_tag,
        )


if __name__ == "__main__":
    app()
