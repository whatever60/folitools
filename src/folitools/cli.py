from typing import Annotated, Literal
from pathlib import Path
import subprocess
import sys
import os

from cyclopts import App, Parameter

from .get_matrix import read_counts
from .primer_info import get_read_stats
from .utils import expand_path_to_list

app = App(help="Foli Tools CLI")
FASTQ_EXTENSIONS = ["fq", "fastq", "fq.gz", "fastq.gz"]
BAM_EXTENSIONS = ["bam", "sam", "sam.gz", "sam.bz2"]


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
        list[str],
        Parameter(
            consume_multiple=True,
            help="File path, glob pattern, or list of file paths for R1 FASTQ files",
        ),
    ],
    output_dir: Annotated[str, Parameter(help="Output directory for trimmed files")],
    cores: Annotated[int, Parameter(help="Number of cores to use")] = 8,
    skip: Annotated[int, Parameter(help="Number of samples to skip")] = 0,
    delete: Annotated[
        bool, Parameter(help="Delete input files after processing")
    ] = False,
) -> None:
    """
    Run fastp preprocessing.

    Args:
        input_: File path, glob pattern, or list of file paths for R1 FASTQ files.
        output_dir: Output directory for trimmed files.
        cores: Number of CPU cores to allocate for fastp.
    """
    # Expand patterns to actual file paths
    file_paths = expand_path_to_list(input_, suffix=FASTQ_EXTENSIONS)
    # Convert list to space-separated string for shell script
    input_patterns = " ".join(file_paths)
    run(
        f"{os.path.dirname(__file__)}/scripts/foli_01_fastp.sh",
        (input_patterns, output_dir, cores, skip, delete),
    )


@app.command(help="Run cutadapt demultiplexing")
def assign_probes(
    *,
    input_: Annotated[
        list[str],
        Parameter(
            consume_multiple=True,
            help="File path, glob pattern, or list of file paths for R1 FASTQ files",
        ),
    ],
    output_dir: Annotated[str, Parameter(help="Output directory for UMI-tagged files")],
    i5: Annotated[Path, Parameter(help="Path to i5 adapter FASTA file")],
    i7: Annotated[Path, Parameter(help="Path to i7 adapter FASTA file")],
    cores: Annotated[int, Parameter(help="Number of cores to use")] = 8,
    skip: Annotated[int, Parameter(help="Number of samples to skip")] = 0,
    delete: Annotated[
        bool, Parameter(help="Delete input files after processing")
    ] = False,
):
    """Run the cutadapt step of the pipeline."""
    # Expand patterns to actual file paths
    file_paths = expand_path_to_list(input_, suffix=FASTQ_EXTENSIONS)
    # Convert list to space-separated string for shell script
    input_patterns = " ".join(file_paths)
    run(
        "foli_02_cutadapt.sh",
        (
            input_patterns,
            output_dir,
            str(i5),
            str(i7),
            str(cores),
            str(skip),
            str(delete),
        ),
    )


@app.command(help="Run mapping step")
def map_(
    *,
    input_: Annotated[
        list[str],
        Parameter(
            consume_multiple=True,
            help="File path, glob pattern, or list of file paths for R1 FASTQ or BAM/SAM files.",
        ),
    ],
    output_bam: Annotated[str, Parameter(help="Output directory for BAM files")],
    output_star: Annotated[
        str | None,
        Parameter(
            help="Output directory for STAR alignment files (required if any FASTQ inputs)"
        ),
    ] = None,
    star_index: Annotated[
        Path | None,
        Parameter(
            "--star-index",
            help="Path to the STAR genome index directory (required if any FASTQ inputs).",
        ),
    ] = None,
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
    strand: Annotated[
        Literal[0, 1, 2],
        Parameter(
            "--strand",
            help="Strand specificity for featureCounts (0=unstranded, 1=stranded, 2=reversely stranded).",
        ),
    ] = 0,
    allow_overlap: Annotated[
        bool,
        Parameter(
            "--allow-overlap",
            help="Allow reads to be assigned to overlapping features. When enabled, passes -O and --fraction to featureCounts.",
        ),
    ] = False,
    allow_multimapping: Annotated[
        bool,
        Parameter(
            "--allow-multimapping",
            help="Allow multi-mapping reads to be counted. When enabled, passes -M and --fraction to featureCounts.",
        ),
    ] = False,
    skip: Annotated[int, Parameter(help="Number of samples to skip")] = 0,
    delete: Annotated[
        bool, Parameter(help="Delete input files after processing")
    ] = False,
):
    """Run the mapping step of the pipeline. Supports both FASTQ (.fq/.fastq/.fq.gz/.fastq.gz) and BAM/SAM inputs. For FASTQ files, includes read filtering to remove short reads (primer dimers) and STAR alignment. For BAM/SAM files, skips alignment and goes directly to featureCounts."""
    # Expand patterns to actual file paths
    file_paths = expand_path_to_list(input_, suffix=FASTQ_EXTENSIONS + BAM_EXTENSIONS)
    # Convert list to space-separated string for shell script
    input_patterns = " ".join(file_paths)
    run(
        "foli_03_map.sh",
        (
            input_patterns,
            output_bam,
            output_star or "",
            str(star_index) if star_index else "",
            str(gtf),
            str(cores),
            str(skip),
            str(delete),
            str(strand),
            str(allow_overlap),
            str(allow_multimapping),
        ),
    )


@app.command(help="Count UMI with umi_tools")
def count(
    *,
    input_: Annotated[
        list[str],
        Parameter(
            consume_multiple=True,
            help="File path, glob pattern, or list of file paths for BAM files",
        ),
    ],
    output_dir: Annotated[str, Parameter(help="Output directory for count results")],
    cores: Annotated[int, Parameter(help="Number of cores to use")] = 8,
    skip: Annotated[int, Parameter(help="Number of samples to skip")] = 0,
):
    """Run the counting step of the pipeline."""
    # Expand patterns to actual file paths
    file_paths = expand_path_to_list(input_, suffix=BAM_EXTENSIONS)
    # Convert list to space-separated string for shell script
    input_patterns = " ".join(file_paths)
    run("foli_04_count.sh", (input_patterns, output_dir, str(cores), str(skip)))


@app.command(help="Generate count matrix")
def get_count_mtx(
    *,
    input_: list[str],
    output: Annotated[str, Parameter(help="Output file path (default: stdout)")],
    gtf: str | None = None,
) -> None:
    """
    Args:
        input_: A file path, glob pattern, or list of file paths.
        gtf: Optional path to a GTF file for gene_id→gene_symbol mapping.
    """
    df = read_counts(input_, gtf)
    if any(output.endswith(ext) for ext in [".tsv", ".txt", ".tsv.gz", ".txt.gz"]):
        sep = "\t"
    elif any(output.endswith(ext) for ext in [".csv", ".csv.gz"]):
        sep = ","
    else:
        raise ValueError(
            f"Output file must end with .tsv, .txt, .tsv.gz, or .txt.gz: {output!r}"
        )
    df.to_csv(output, sep=sep, index_label="gene", header=True)


@app.command(help="Get read statistics from FASTQ files after primer assignment")
def get_read_stats_(
    *,
    input_: Annotated[
        list[str],
        Parameter(
            consume_multiple=True,
            help="File path, glob pattern, or list of file paths for R1 FASTQ files",
        ),
    ],
    output_dir: str,
    cores: int = 1,
    overwrite: bool = False,
    skip: Annotated[int, Parameter(help="Number of samples to skip")] = 0,
) -> None:
    """
    Process FASTQ files to extract read statistics and write to parquet files.

    Args:
        input_: A file path, glob pattern, or list of file paths for R1 FASTQ files.
        output_dir: Directory to save the output parquet files.
        cores: Number of CPU cores to use for parallel processing.
        overwrite: If True, overwrite existing parquet files.
    """
    get_read_stats(
        read1_pattern=input_,
        output_dir=output_dir,
        cores=cores,
        overwrite=overwrite,
        skip=skip,
    )


if __name__ == "__main__":
    app()
