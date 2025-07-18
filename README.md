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

All the output directories (`./fastp/`, `./rest/`, `./featurecounts/`, etc.) are automatically created in the current working directory if they don't exist.

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

Note that input FASTQ files are assumed to be paired-end. Read 1 pattern is provided to the script by users and files of read 2 will be automatically derived from read 1. This is also true for the following steps.

Optionally, you can restrict processing to a subset of samples by providing a custom glob pattern.

For example:
```bash
foli fastp 'sample-A*_R1_*.fastq.gz'
```

### Step 2. Probe assignment

```bash
foli assign-probes
```

This command performs primer or adapter trimming using `cutadapt`. It processes reads from the `./fastp/` directory and writes output FASTQs to `./rest/` and `./rest_all/`. Primer or barcode sequences should be provided as FASTA files in `./data/` by default.

Optionally, you can:
- Specify a custom glob pattern to restrict which files are processed (e.g. `'sample_A*_1.fq.gz'`)
- Provide a different adapter directory (e.g. `./barcodes`)

For example:
```bash
foli cutadapt --pattern 'sample-A*_1.fq.gz' --adapter-dir ./barcodes --threads 8
```

### Step 3. Mapping and Feature Counting

```bash
foli map --star-index STAR_INDEX_PATH --gtf GTF_PATH
```

This command aligns reads to a genome using `STAR` and assigns transcripts using `featureCounts`. It expects paired-end FASTQs in `./rest/`, and writes `STAR` output to `./star/` and `featureCounts` output (sorted BAM files and count tables) to `./featurecounts/`.

Optionally, you can provide a glob pattern to filter input files.

For example:
```bash
foli map --pattern 'sample-A*_1.fq.gz' --star-index STAR_INDEX_PATH --gtf GTF_PATH
```

### Step 4. Gene Counting

```bash
foli count
```

This command runs `umi_tools` on BAM files in `./featurecounts/` and outputs read count tables to `./counts/`.

Optionally, you can provide a glob pattern to restrict input BAMs.

For example:
```bash
foli count --pattern 'sample-A*.bam'
```

## Output Structure

Each stage writes outputs to stage-specific subdirectories:

| Stage      | Output Directories                    | Key Files |
|------------|---------------------------------------|-----------|
| fastp      | `./fastp/`, `./fastp_fastqc/`        | Trimmed FASTQ files, QC reports |
| cutadapt   | `./rest_all/`, `./rest/`             | Adapter-trimmed reads, UMI-tagged |
| map        | `./star/`, `./featurecounts/`        | Alignments, sorted BAM files |
| count      | `./counts/`                          | UMI count matrices, grouped BAMs |
