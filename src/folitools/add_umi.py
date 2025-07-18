import sys
import gzip
from typing import Generator, TextIO, Annotated

from cyclopts import Parameter, run


def read_interleaved_fastq(
    fastq_stream,
) -> Generator[tuple[str, str, str, str, str, str], None, None]:
    while True:
        header_1 = fastq_stream.readline().strip()
        if not header_1:
            break
        seq_1 = fastq_stream.readline().strip()
        fastq_stream.readline()  # +
        qual_1 = fastq_stream.readline().strip()

        header_2 = fastq_stream.readline().strip()
        seq_2 = fastq_stream.readline().strip()
        fastq_stream.readline()  # +
        qual_2 = fastq_stream.readline().strip()

        yield header_1, seq_1, qual_1, header_2, seq_2, qual_2


def add_umi(
    fastq_stream: TextIO,
    output_fastq_1: str,  # assumed to be gzipped
    output_fastq_2: str,  # assumed to be gzipped
    sep: str = "_",
):
    with (
        gzip.open(output_fastq_1, "wt") as out_fastq_1,
        gzip.open(output_fastq_2, "wt") as out_fastq_2,
    ):
        for (
            header_1,
            seq_1,
            qual_1,
            header_2,
            seq_2,
            qual_2,
        ) in read_interleaved_fastq(fastq_stream):
            read_id_1, adapter_1, *match_sequence_1 = header_1.split(" ")
            read_id_2, adapter_2, *match_sequence_2 = header_2.split(" ")

            if read_id_1 != read_id_2:
                raise ValueError("FASTQ records are not interleaved")

            if match_sequence_1:
                match_sequence_1 = match_sequence_1[0]
                umi_1 = seq_1.split(match_sequence_1)[0]
            else:
                match_sequence_1 = ""
                umi_1 = ""
            if match_sequence_2:
                match_sequence_2 = match_sequence_2[0]
                umi_2 = seq_2.split(match_sequence_2)[0]
            else:
                match_sequence_2 = ""
                umi_2 = ""

            new_id_1 = sep.join([read_id_1, umi_1, umi_2])
            new_id_2 = sep.join([read_id_2, umi_1, umi_2])
            # new_comment_1 = adapter_1
            # new_comment_2 = adapter_2
            new_comment = f"{adapter_1}+{adapter_2}"

            new_seq_1 = seq_1[len(umi_1) :]
            new_seq_2 = seq_2[len(umi_2) :]
            new_qual_1 = qual_1[len(umi_1) :]
            new_qual_2 = qual_2[len(umi_2) :]

            print(
                f"{new_id_1} {new_comment}",
                new_seq_1,
                "+",
                new_qual_1,
                sep="\n",
                file=out_fastq_1,
            )
            print(
                f"{new_id_2} {new_comment}",
                new_seq_2,
                "+",
                new_qual_2,
                sep="\n",
                file=out_fastq_2,
            )


def main(
    o1: Annotated[
        str,
        Parameter(help="Output FASTQ file for read 1 (gzipped)"),
    ],
    o2: Annotated[
        str,
        Parameter(help="Output FASTQ file for read 2 (gzipped)"),
    ],
):
    """
    Add UMI to FASTQ records from stdin and write to output files.

    Args:
        o1: Output FASTQ file for read 1 (gzipped).
        o2: Output FASTQ file for read 2 (gzipped).
        sep: Separator between read ID and UMIs.
    """
    add_umi(fastq_stream=sys.stdin, output_fastq_1=o1, output_fastq_2=o2)


if __name__ == "__main__":
    run(main)
