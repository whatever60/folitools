from typing import Annotated
from pathlib import Path
import subprocess
import sys
import os

from cyclopts import App, Parameter


app = App(help="Foli Tools CLI")

# Get the directory where the scripts are located, relative to this CLI script.


def run(script_name: str, args: tuple) -> None:
    script_path = Path(__file__).parent / "scripts" / script_name
    if not script_path.exists():
        sys.exit(f"ERROR: script not found: {script_path}")
    try:
        subprocess.run(["bash", str(script_path), *list(map(str, args))], check=True)
    except subprocess.CalledProcessError as e:
        sys.exit(
            f"ERROR: Script failed with exit code {e.returncode}\n{e.stderr or ''}"
        )


@app.command(help="Run fastp preprocessing")
def qc(
    *,
    input_: Annotated[
        str, Parameter(help="pattern for input FASTQ files")
    ] = "*_R1_*.fastq.gz",
    cores: Annotated[int, Parameter(help="Number of cores to use")] = 8,
) -> None:
    """
    Run fastp preprocessing.

    Args:
        cores: Number of CPU cores to allocate for fastp.
    """
    # pass the core count as the only argument to your shell script
    run(f"{os.path.dirname(__file__)}/scripts/foli_01_fastp.sh", (input_, cores))


@app.command(help="Run cutadapt demultiplexing")
def assign_probes(
    *,
    input_: Annotated[
        str, Parameter(help="pattern for R1 FASTQ files from fastp output")
    ] = "*_1.fq.gz",
    adapter_dir: Annotated[
        Path, Parameter(help="Directory containing adapter FASTA files")
    ] = Path("./data"),
    cores: Annotated[int, Parameter(help="Number of cores to use")] = 8,
):
    """Run the cutadapt step of the pipeline."""
    run("foli_02_cutadapt.sh", (input_, adapter_dir, str(cores)))


@app.command(help="Run mapping step")
def map_(
    *,
    input_: Annotated[
        str,
        Parameter(
            "--input_",
            help="Optional pattern for R1 FASTQ files (default: '*_1.fq.gz').",
        ),
    ] = "*_1.fq.gz",
    star_index: Annotated[
        Path,
        Parameter(
            "--star-index",
            help="Path to the STAR genome index directory.",
        ),
    ],
    gtf: Annotated[
        Path,
        Parameter(
            "--gtf",
            help="Path to the GTF annotation file for featureCounts.",
        ),
    ],
    cores: Annotated[
        int,
        Parameter(
            "--cores",
            help="Total number of cores to allocate for the pipeline.",
        ),
    ] = 8,
):
    """Run the mapping step of the pipeline. Includes read filtering to remove short reads (primer dimers)."""
    run(
        "foli_03_map.sh",
        (
            "--cores",
            str(cores),
            "--star-index",
            str(star_index),
            "--gtf",
            str(gtf),
            "--input_",
            input_,
        ),
    )


@app.command(help="Run read counting")
def count(
    *,
    input_: Annotated[str, Parameter(help="pattern for BAM files")] = "*.bam",
    cores: Annotated[int, Parameter(help="Number of cores to use")] = 8,
):
    """Run the counting step of the pipeline."""
    run("foli_04_count.sh", (input_, str(cores)))


if __name__ == "__main__":
    app()
