# folitools

**folitools** is a modular CLI toolkit and python library for Foli-seq data processing.

Foli-seq is a amplicon-based high-throughput sequencing method for profiling host gene expression from feces. Foli-seq is developed to probe amplicons of 320-380bp and sequence on PE150 Illumina sequencers.

---

## Installation

Install from PyPI:

```bash
pip install folitools
```

The dependencies from conda are also required:

```bash
conda install -c bioconda -c conda-forge \
  fastp cutadapt samtools bwa-mem2 star seqkit fastqc subread sambamba pigz
```

## Usage

Each stage of the pipeline is exposed as a command via `foli <subcommand>`.

**A note on input paths**: All commands that accept input paths through the `--input` parameter allow files locations of absolute paths, relative paths, or glob patterns. Multiple paths/patterns are also allowed.

**A note on paired-end files**: Foli-seq uses paired-end sequencing. When inputing FASTQ files, only R1 path is given to the command and the corresponding R2 files are automatically derived by replacing the R1 pattern with R2 pattern (e.g., `_R1_` → `_R2_`, `_1` → `_2`). Both R1 and R2 files are required for the pipeline to work.

### Expected Inputs

- Raw paired-end FASTQ files (following the naming pattern `*_R1_*.fastq.gz` and `*_R2_*.fastq.gz`)
- Adapter FASTA files
- STAR genome index directory
- GTF annotation file

Output directories are automatically created if they don't exist.

### Step 1. Preprocessing

```bash
foli qc --input "/path/to/fastq/*_R1_*.fastq.gz" --output-dir "/path/to/trimmed_reads"
```

This command performs read trimming using `fastp` and runs quality check on output files with `fastqc` and `seqkit`.

You can specify input FASTQ files from any location using absolute paths, relative paths, or glob patterns. The paired-end R2 files are automatically derived from R1 files by replacing `_R1_` with `_R2_`.

Output files:
- Trimmed reads in the specified output directory
- FastQC reports for trimmed reads
- Summary statistics for trimmed reads

### Step 2. Probe assignment

```bash
foli assign-probes --input "/path/to/trimmed_reads/*_1.fq.gz" --output-dir "/path/to/umi_tagged" --i5 /path/to/i5_short.fasta --i7 /path/to/i7_short.fasta
```

This command performs gene probe primer assignment and trimming using `cutadapt`. This step extracts UMI sequences from adapter matches and embeds them in read names for downstream processing. Quality statistics are also generated for the output files using `seqkit`.

After probe assignment, you can get read statistics:

```bash
foli get-read-stats --input "/path/to/umi_tagged/*_1.fq.gz" --output-dir "/path/to/read_stats"
```

This analyzes the UMI-tagged reads to generate detailed statistics about primer assignment and read characteristics.

Output files:
- UMI-tagged, adapter-trimmed reads
- Summary statistics for UMI-tagged reads
- Read statistics in parquet format

### Step 3. Mapping and Feature Counting

```bash
foli map --input "/path/to/umi_tagged/*_1.fq.gz" --output-bam "/path/to/mapped_reads" --output-star "/path/to/star" --star-index /path/to/star_index --gtf /path/to/annotation.gtf
```

This command performs streamlined alignment and feature counting using `STAR` and `featureCounts`. This step filters out short reads (< 60bp, considered primer dimers) using `cutadapt`, aligns filtered reads to the genome using `STAR`, and assigns reads to genomic features using `featureCounts`.

UMI sequences and cell barcodes are added to BAM records as tags (`US` for raw UMI, `UC` for filtered UMI, `CB` for cell barcode) for the next counting step.

Output files:
- STAR alignment files in the specified STAR output directory
- Sorted, tagged BAM files in the specified BAM output directory
- `featureCounts` tables and logs

### Step 4. UMI-based Gene Counting

```bash
foli count --input "/path/to/mapped_reads/*.sorted.bam" --output-dir "/path/to/count_results"
```

This command processes BAM files using `umi_tools group` to generate UMI-deduplicated count data. This step:

1. Groups reads by UMI, cell barcode, and gene assignment
2. Handles paired-end reads, unmapped reads, and chimeric pairs
3. Outputs detailed grouping information for downstream analysis

Input BAM files should contain UMI tags (`UC`), cell tags (`CB`), and gene assignment tags (`XT`) from the previous step.

After counting, generate the final count matrix:

```bash
foli get-count-mtx --input "/path/to/count_results/*.group.tsv.gz" --output "/path/to/final_counts.tsv" --gtf /path/to/annotation.gtf
```

You can also use the Python function directly:

```python
from folitools.get_matrix import read_counts
df = read_counts("/path/to/input_files", "/path/to/annotation.gtf")
df.to_csv("/path/to/output.tsv", sep="\t", index_label="gene")
```

Output files:
- UMI grouping tables (`.group.tsv.gz`)
- Processing logs for each sample
- Final count matrix in TSV format


## TODO

- Add QC report code that summarizes all metrics into a table (need refactoring)
- Add probe set design code
- Primer set evaluation (fasta recover and primer dimer)
