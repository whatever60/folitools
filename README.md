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

⚠️ The current working directory (`pwd`) is treated as the working environment. All outputs will be written to subdirectories (e.g. `./fastp/`, `./featurecounts/`, `./counts/`), and required input files must be present or referenced relative to this directory.

### Expected Inputs

- `./fastq/`: directory containing raw paired-end FASTQ files
  - Files should follow the naming pattern `*_R1_*.fastq.gz` and `*_R2_*.fastq.gz`
- `./data/`: directory containing adapter FASTA files (named as `i5_short.fasta` and `i7_short.fasta`)
- STAR genome index directory
- GTF annotation file

All the output directories (`./fastp/`, `./rest_all/`, `./featurecounts/`, etc.) are automatically created in the current working directory if they don't exist.

### Step 1. Preprocessing

```bash
foli qc
```

This command performs quality control and read trimming using `fastqc`, `seqkit`, and `fastp`.

Input FASTQ files are expected in `./fastq/`, with paired-end naming like `*_R1_*.fastq.gz` and corresponding R2 files.
Outputs include:
- Trimmed reads in `./fastp/`
- FastQC reports in `./fastq_fastqc/` and `./fastp_fastqc/`
- Summary statistics in `fastq.stats` and `fastp.stats`

Note that input FASTQ files are assumed to be paired-end. Read 1 pattern is provided to the script by users and files of read 2 will be automatically derived from read 1. This is also the case for the following steps.

Optionally, you can restrict processing to a subset of samples by providing a custom glob pattern.

For example:
```bash
foli qc --input 'sample-A*_R1_*.fastq.gz' --cores 8
```

### Step 2. Probe assignment

```bash
foli assign-probes
```

This command performs gene probe primer assignment and trimming using `cutadapt`. It processes reads from the `./fastp/` directory and writes output FASTQs to `./rest_all/`.

Probe primer sequences should be provided as FASTA files in `./data/` by default.

This step extracts UMI sequences from adapter matches and embeds them in read names for downstream processing. A summary statistics of FASTQ files in `./rest_all` is also generated at `./rest_all.stats`.

Optionally, you can:
- Specify a custom glob pattern to restrict which files are processed (e.g. `'sample_A*_1.fq.gz'`)
- Provide a different adapter directory (e.g. `./barcodes`)

For example:
```bash
foli assign-probes --input 'sample-A*_1.fq.gz' --adapter-dir ./barcodes --cores 8
```

### Step 3. Mapping and Feature Counting

```bash
foli map --star-index STAR_INDEX_PATH --gtf GTF_PATH
```

This command performs a streamlined alignment and feature counting using `STAR` and `featureCounts`.

This step filters out short reads (< 60bp, considered primer dimers) using `cutadapt`, aligns filtered reads to the genome using `STAR`, and assigns reads to genomic features using `featureCounts`.

UMI sequences and cell barcodes are added to BAM records as tags (`US` for raw UMI, `UC` for filtered UMI, `CB` for cell barcode) for the next counting step.

Input files are expected in `./rest_all/` with paired-end naming pattern.

Outputs include:
- STAR alignment files in `./star/`
- Sorted, tagged BAM files in `./featurecounts/`
- `featureCounts` tables and logs

Optionally, you can provide a glob pattern to filter input files.

For example:
```bash
foli map --input 'sample-A*_1.fq.gz' --star-index STAR_INDEX_PATH --gtf GTF_PATH --cores 8
```

### Step 4. UMI-based Gene Counting

```bash
foli count
```

This command processes BAM files from `./featurecounts/` using `umi_tools group` to generate UMI-deduplicated count data. This step:

1. Groups reads by UMI, cell barcode, and gene assignment
2. Handles paired-end reads, unmapped reads, and chimeric pairs
3. Outputs detailed grouping information for downstream analysis

Input BAM files are expected to contain UMI tags (`UC`), cell tags (`CB`), and gene assignment tags (`XT`) from the previous step.

Outputs include:
- UMI grouping tables (`.group.tsv.gz`) in `./counts/`
- Processing logs for each sample

Optionally, you can provide a glob pattern to restrict input BAMs.

For example:
```bash
foli count --input 'sample-A*.bam' --cores 8
```

## Output Structure

Each stage writes outputs to stage-specific subdirectories:

| Stage         | Output Directories                  | Key Files                           |
|---------------|-------------------------------------|-------------------------------------|
| qc            | `./fastp/`, `./fastp_fastqc/`       | Trimmed FASTQ files, QC reports     |
| assign-probes | `./rest_all/`                       | UMI-tagged, adapter-trimmed reads   |
| map           | `./star/`, `./featurecounts/`       | Alignments, sorted tagged BAM files |
| count         | `./counts/`                         | UMI grouping tables, logs           |

## TODO

- Add QC report code
- Add count matrix generation code and basic analysis