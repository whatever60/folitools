import sys
import io
from typing import Annotated

import pysam
from cyclopts import App, Parameter


app = App()


def add_tags_wo_fastq(
    bam_input: str,
    bam_output: str,
    umi_tag_s_name: str = "US",  # single-read UMI
    umi_tag_p_name: str = "UP",  # concatenated UMI from paired reads
    cell_tag_name: str = "CB",
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
        # Peek at stdin to determine if it's BAM or SAM
        initial_bytes = sys.stdin.buffer.read(4)
        is_bam = initial_bytes.startswith(b"BAM\x01")
        # Try to seek back if possible
        if sys.stdin.buffer.seekable():
            sys.stdin.buffer.seek(-len(initial_bytes), 1)
            input_ = sys.stdin.buffer
        else:
            # Create a file-like object that prepends the already-read bytes
            input_ = io.BytesIO(initial_bytes + sys.stdin.buffer.read())
    else:  # An actual file path
        with open(bam_input, "rb") as f:
            is_bam = f.read(4).startswith(b"BAM\x01")
        input_ = bam_input
    mode_in = "rb" if is_bam else "r"

    mode_out = "w"
    if bam_output == "-":
        mode_out = "wb0"  # same as wbu, i.e., uncompressed BAM
    elif bam_output.endswith("bam"):
        mode_out = "wb"
    elif bam_output.endswith("sam"):
        mode_out = "w"
    else:
        raise ValueError("Output file must end with .bam or .sam")

    with (
        pysam.AlignmentFile(input_, mode_in) as sam_in,
        pysam.AlignmentFile(bam_output, mode_out, header=sam_in.header) as bam_out,
    ):
        for read in sam_in:
            query_name = read.query_name
            assert query_name is not None, "Missing query name in read"
            read_name, seq1, seq2 = query_name.split("_")
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
            print(query_name, file=sys.stderr)
            bam_out.write(read)


@app.command(help="Add UMI and cell tags to BAM file from FASTQ sequences.")
def main(
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
) -> None:
    """
    Dispatch to add_tags or add_tags_wo_fastq based on presence of --index_fastq.

    All existing comments and logic are preserved.
    """
    add_tags_wo_fastq(
        bam_input="-",
        bam_output=output,
        cell_tag_name=cell_tag_name,
        cell_tag=cell_tag,
    )


if __name__ == "__main__":
    app()
