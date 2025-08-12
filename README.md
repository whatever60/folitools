# folitools

**folitools** is a lightweight, modular RNA-seq preprocessing toolkit built from shell scripts and wrapped with a Python CLI.

---

## Installation

Install from PyPI:

```bash
pip install folitools
```

Then manually install the required bioinformatics tools via conda:

```bash
conda install -c bioconda -c conda-forge \
  fastp cutadapt samtools bwa-mem2 star seqkit fastqc subread sambamba pigz
```

## Usage

Each stage of the pipeline is exposed as a command via `foli <subcommand>`, which wraps the corresponding shell script. The CLI uses explicit arguments with proper help documentation.

All commands accept flexible input paths through the `--input` parameter, allowing you to specify files from any location using absolute paths, relative paths, or glob patterns. You can also provide multiple patterns or file paths by specifying the `--input` parameter multiple times. Additionally, all commands now support custom output directories through the `--output-dir` parameter for better organization.

**Note on paired-end files**: All input specifications should refer to R1 files. The corresponding R2 files are automatically derived by replacing the R1 pattern with R2 pattern (e.g., `_R1_` → `_R2_`, `_1` → `_2`).

### Expected Inputs

- Raw paired-end FASTQ files (following the naming pattern `*_R1_*.fastq.gz` and `*_R2_*.fastq.gz`)
- Adapter FASTA files (named as `i5_short.fasta` and `i7_short.fasta`) in a directory specified by `--probe-dir`
- STAR genome index directory
- GTF annotation file

Output directories are automatically created if they don't exist.

### Step 1. Preprocessing

```bash
foli qc --input "/path/to/fastq/*_R1_*.fastq.gz" --output-dir "./trimmed_reads"
```

This command performs quality control and read trimming using `fastqc`, `seqkit`, and `fastp`.

You can specify input FASTQ files from any location using absolute paths, relative paths, or glob patterns. The paired-end R2 files are automatically derived from R1 files by replacing `_R1_` with `_R2_`.

Outputs include:
- Trimmed reads in the specified output directory (default: `./fastp/`)
- FastQC reports in `{output_dir}_fastqc/` and `./fastq_fastqc/`
- Summary statistics files

You can specify multiple patterns or individual files:

```bash
foli qc --input "/data/batch1/*_R1_*.fastq.gz" "/data/batch2/sample_*_R1_*.fastq.gz" --output-dir "./my_trimmed" --cores 8
```

### Step 2. Probe assignment

```bash
foli assign-probes --input "./trimmed_reads/*_1.fq.gz" --output-dir "./umi_tagged" --probe-dir ./data
```

This command performs gene probe primer assignment and trimming using `cutadapt`. You can specify input files from any location, and both the output and probe directories can be customized.

This step extracts UMI sequences from adapter matches and embeds them in read names for downstream processing. Output FASTQs are written to the specified directory, and summary statistics are generated.

You can specify multiple input patterns and customize directories:

```bash
foli assign-probes --input "/data/fastp/sample-A*_1.fq.gz" "/data/fastp/sample-B*_1.fq.gz" --output-dir "./my_umi_files" --probe-dir ./custom_probes --cores 8
```

### Step 3. Mapping and Feature Counting

```bash
foli map --input "./umi_tagged/*_1.fq.gz" --output-dir "./mapped_reads" --star-index STAR_INDEX_PATH --gtf GTF_PATH
```

This command performs streamlined alignment and feature counting using `STAR` and `featureCounts`. You can specify input files from any location and customize the output directory.

This step filters out short reads (< 60bp, considered primer dimers) using `cutadapt`, aligns filtered reads to the genome using `STAR`, and assigns reads to genomic features using `featureCounts`.

UMI sequences and cell barcodes are added to BAM records as tags (`US` for raw UMI, `UC` for filtered UMI, `CB` for cell barcode) for the next counting step.

Outputs include:
- STAR alignment files in the parent directory of output_dir (e.g., `./star/` if output_dir is `./mapped_reads`)
- Sorted, tagged BAM files in the specified output directory
- `featureCounts` tables and logs

You can specify multiple input patterns:

```bash
foli map --input "/data/umi_files/sample-A*_1.fq.gz" --output-dir "./my_mapped" --star-index STAR_INDEX_PATH --gtf GTF_PATH --cores 8
```

### Step 4. UMI-based Gene Counting

```bash
foli count --input "./mapped_reads/*.sorted.bam" --output-dir "./count_results"
```

This command processes BAM files using `umi_tools group` to generate UMI-deduplicated count data. You can specify BAM files from any location and customize the output directory.

This step:

1. Groups reads by UMI, cell barcode, and gene assignment
2. Handles paired-end reads, unmapped reads, and chimeric pairs
3. Outputs detailed grouping information for downstream analysis

Input BAM files should contain UMI tags (`UC`), cell tags (`CB`), and gene assignment tags (`XT`) from the previous step.

Outputs include:
- UMI grouping tables (`.group.tsv.gz`) in the specified directory
- Processing logs for each sample

You can specify multiple input patterns:

```bash
foli count --input "/data/featurecounts/sample-A*.bam" "/data/featurecounts/sample-B*.bam" --output-dir "./my_counts" --cores 8
```

### Complete Pipeline Example

Here's a complete example showing how to run the entire pipeline with flexible input and output paths:

```bash
# Set environment variables
export STAR_INDEX="/path/to/star/index"
export GTF_FILE="/path/to/annotation.gtf"

# Step 1: Quality control and preprocessing
foli qc --input "/data/raw_fastq/*_R1_001.fastq.gz" --output-dir "./trimmed_reads" --cores 16

# Step 2: Probe assignment
foli assign-probes --input "./trimmed_reads/*_1.fq.gz" --output-dir "./umi_tagged" --probe-dir /path/to/probes --cores 16

# Step 3: Get read statistics (optional)
foli get-read-stats --input "./umi_tagged/*_1.fq.gz" --output-dir "./read_stats" --cores 16

# Step 4: Mapping and feature counting
foli map --input "./umi_tagged/*_1.fq.gz" --output-dir "./mapped_reads" --star-index "$STAR_INDEX" --gtf "$GTF_FILE" --cores 16

# Step 5: UMI-based counting
foli count --input "./mapped_reads/*.sorted.bam" --output-dir "./count_results" --cores 16

# Step 6: Generate count matrix
foli get-count-mtx --input "./count_results/*.group.tsv.gz" --output "./final_counts.tsv" --gtf "$GTF_FILE"
```

### Multiple Input Examples

You can specify multiple input patterns or files:

```bash
# Multiple glob patterns with custom output directories
foli qc --input "/batch1/*_R1_*.fastq.gz" "/batch2/*_R1_*.fastq.gz" --output-dir "./all_trimmed" --cores 16

# Mix of patterns and specific files
foli map --input "./sample_A*_1.fq.gz" "./sample_B_special_1.fq.gz" --output-dir "./custom_mapped" --star-index "$STAR_INDEX" --gtf "$GTF_FILE"

# Using absolute paths from different locations
foli count --input "/project1/bams/*.bam" "/project2/bams/*.bam" --output-dir "/results/combined_counts" --cores 8
```

## Output Structure

Each stage writes outputs to customizable directories (defaults shown):

| Stage         | Default Output Directories         | Key Files                           |
|---------------|-------------------------------------|-------------------------------------|
| qc            | `./fastp/`, `./fastp_fastqc/`       | Trimmed FASTQ files, QC reports     |
| assign-probes | `./rest_all/`                       | UMI-tagged, adapter-trimmed reads   |
| map           | `./featurecounts/`, `./star/`       | Alignments, sorted tagged BAM files |
| count         | `./counts/`                         | UMI grouping tables, logs           |

All output directories can be customized using the `--output-dir` parameter.

## TODO

- Add QC report code (need refactoring)
- Add probe set design code
- Rarefying curves
- Primer set evaluation
- Move seqkit, fastQC and STAR index loading to their own modules?