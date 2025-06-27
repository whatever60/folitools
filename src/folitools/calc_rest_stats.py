import glob
import gzip
import os
from pathlib import Path
from typing import Annotated

import polars as pl
from Bio import SeqIO
from joblib import Parallel, delayed
from cyclopts import App, Parameter
from tqdm.auto import tqdm


BATCH_SIZE = 1_000_000

app = App(help_format="plaintext")


def process_sample(read1_path: str, output_dir: str) -> None:
    """
    Process a single pair of FASTQ files (R1 and R2), extract UMI and primer info,
    and write batched intermediate results to temp files, combining at the end.

    Args:
        read1_path: Path to the R1 FASTQ file (.fq.gz).
        output_dir: Directory to write the final Parquet output.
    """
    sample = Path(read1_path).name.split("_")[0]
    out_path = Path(output_dir) / f"{sample}.parquet"

    if os.path.exists(out_path):
        print(f"Output file {out_path} already exists. Skipping sample {sample}.")
        return

    tmp_dir = Path(output_dir) / "tmp"
    os.makedirs(tmp_dir, exist_ok=True)

    read2_path = read1_path.replace("_1.fq.gz", "_2.fq.gz")
    buffer = []
    batch_idx = 0
    temp_paths = []

    with (
        gzip.open(read1_path, "rt") as r1_handle,
        gzip.open(read2_path, "rt") as r2_handle,
    ):
        for r1, r2 in tqdm(
            zip(SeqIO.parse(r1_handle, "fastq"), SeqIO.parse(r2_handle, "fastq")),
            desc=f"Reading {sample}",
        ):
            if r1.id != r2.id:
                raise ValueError(
                    f"Read ID mismatch in sample {sample}: {r1.id} vs {r2.id}"
                )
            primer_fwd, primer_rev = r1.description.split(" ")[-1].split("+")
            read_id, umi5, umi3 = r1.id.split("_")
            buffer.append(
                {
                    "read_id": read_id,
                    "primer_fwd": primer_fwd,
                    "primer_rev": primer_rev,
                    "r1_length": len(r1.seq),
                    "r2_length": len(r2.seq),
                    "umi5": umi5,
                    "umi3": umi3,
                }
            )

            if len(buffer) >= BATCH_SIZE:
                tmp_path = tmp_dir / f"{sample}_batch{batch_idx}.parquet"
                pl.DataFrame(buffer).write_parquet(tmp_path)
                temp_paths.append(tmp_path)
                buffer.clear()
                batch_idx += 1

    if buffer:
        tmp_path = tmp_dir / f"{sample}_batch{batch_idx}.parquet"
        pl.DataFrame(buffer).write_parquet(tmp_path)
        temp_paths.append(tmp_path)

    if not temp_paths:
        print(f"No records for sample {sample}. Skipping.")
        return

    os.makedirs(output_dir, exist_ok=True)
    if len(temp_paths) == 1:
        temp_paths[0].rename(out_path)
    else:
        full_df = pl.concat([pl.read_parquet(p) for p in temp_paths])
        full_df.write_parquet(out_path)
        for p in temp_paths:
            p.unlink()


@app.default
def main(
    input_dir: Annotated[
        str,
        Parameter(
            name=["--input-dir", "-i"],
            help="Directory containing fastq files for read 1.",
        ),
    ],
    output_dir: Annotated[
        str,
        Parameter(
            name=["--output-dir", "-o"], help="Directory to write output .parquet files"
        ),
    ],
    num_cpus: Annotated[
        int,
        Parameter(name="--num-cpus", help="Number of CPUs for parallel processing"),
    ] = 1,
) -> None:
    """
    Extract read stats from paired-end FASTQ files with UMI and primer annotations.
    """
    input_dir_path = os.path.expanduser(input_dir)
    output_dir_path = os.path.expanduser(output_dir)
    read1_files = sorted(glob.glob(os.path.join(input_dir_path, "*_1.fq.gz")))

    Parallel(n_jobs=num_cpus)(
        delayed(process_sample)(read1_path, output_dir_path)
        for read1_path in tqdm(read1_files, desc="Processing samples")
    )


if __name__ == "__main__":
    app()
