import sys
import gzip
import argparse


def read_interleaved_fastq(fastq_stream) -> tuple[str, str, str, str, str, str]:
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
    fastq_stream: sys.stdin,
    output_fastq_1: str,  # assumed to be gzipped
    output_fastq_2: str,  # assumed to be gzipped
    sep="_",
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
            new_comment_1 = adapter_1
            new_comment_2 = adapter_2

            new_seq_1 = seq_1[len(umi_1) :]
            new_seq_2 = seq_2[len(umi_2) :]
            new_qual_1 = qual_1[len(umi_1) :]
            new_qual_2 = qual_2[len(umi_2) :]

            print(
                f"{new_id_1} {new_comment_1}",
                new_seq_1,
                "+",
                new_qual_1,
                sep="\n",
                file=out_fastq_1,
            )
            print(
                f"{new_id_2} {new_comment_2}",
                new_seq_2,
                "+",
                new_qual_2,
                sep="\n",
                file=out_fastq_2,
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Add UMI to FASTQ records from stdin and write to output files."
    )
    parser.add_argument(
        "--o1",
        type=str,
        required=True,
        help="Output FASTQ file for read 1.",
    )
    parser.add_argument(
        "--o2",
        type=str,
        required=True,
        help="Output FASTQ file for read 2.",
    )
    args = parser.parse_args()
    add_umi(
        fastq_stream=sys.stdin,
        output_fastq_1=args.o1,
        output_fastq_2=args.o2,
    )
