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
  fastp cutadapt samtools bwa-mem2 star seqkit fastqc subread pigz
```

## Usage

Each stage of the pipeline is exposed as a command via `foli <stage>`, which wraps the corresponding shell script. Arguments are passed directly to the underlying `.sh` scripts and behave as positional arguments (not keyword-style flags).

⚠️ The current working directory (`pwd`) is treated as the working environment. All outputs will be written to subdirectories (e.g. `./fastp/`, `./star_bam/`, `./counts/`), and required input files must be present or referenced relative to this directory.

### Expected Inputs

- `./fastq/`: directory containing raw paired-end FASTQ files
  - Files should follow the naming pattern `*_R1_*.fastq.gz` and `*_R2_*.fastq.gz`
- `primer.fasta` or similar: barcode or primer file (used in `cutadapt`)
- `transcript.fa`: reference transcriptome FASTA (for index building or annotation)
- `annotation.gtf`: gene annotation file (used in `foli count`)
- STAR genome index directory: passed as an argument to `foli map`

All other intermediate and output directories (`./fastp/`, `./rest/`, `./star/`, etc.) are automatically created if they don’t exist.

### Step 1. Preprocessing

```bash
foli fastp
```

This command performs quality control and read trimming using `fastqc`, `seqkit`, and `fastp`.

Input FASTQ files are expected in `./fastq/`, with paired-end naming like `*_R1_*.fastq.gz` and corresponding R2 files.
Outputs include:
- Trimmed reads in `./fastp/`
- FastQC reports in `./fastq_fastqc/` and `./fastp_fastqc/`
- Summary statistics in `fastq.stats` and `fastp.stats`

You can optionally restrict processing to a subset of samples by providing a custom glob pattern:

```bash
foli fastp 'prefixabc*_R1_*.fastq.gz'
```

### Step 2. Primer Trimming

```bash
foli cutadapt
```

This command performs primer or adapter trimming using `cutadapt`. It processes reads from the `./fastp/` directory and writes output FASTQs to `./rest/` and `./rest_all/`. Primer or barcode sequences should be provided as FASTA files in `./data/` by default.

You can optionally:
- Specify a custom glob pattern to restrict which files are processed (e.g. `'sample*_1.fq.gz'`)
- Provide a different adapter directory (e.g. `./barcodes/`)

Example:
```bash
foli cutadapt 'sample*_1.fq.gz' ./barcodes/
```

### Step 3. Mapping

```bash
foli map STAR_INDEX_PATH
```

This command aligns reads to a genome using STAR. It expects paired-end FASTQs in `./rest/`, and writes sorted BAM files to `./star_bam/`.

The first argument must be the path to an existing STAR index directory.

You can optionally provide a glob pattern to filter input reads:

```bash
foli map 'sample*_1.fq.gz' ./star_index/
```

### Step 4. Gene Counting

```bash
foli count GTF_PATH
```

This command runs `featureCounts` on BAM files in `./star_bam/` and outputs read count tables to `./featurecounts/` and `./counts/`.

The first argument is the path to a GTF annotation file.

To restrict input BAMs, you can provide an optional glob pattern as the first argument and move the GTF path to the second position:

```bash
foli count 'sample*.bam' ./annotation.gtf
```

## Output

Each stage writes outputs to stage-specific subdirectories:

| Stage      | Output Directory            |
|------------|-----------------------------|
| fastp      | `./fastp/`, `./fastp_fastqc/` |
| cutadapt   | `./rest/`, `./rest_all/`     |
| map        | `./star/`, `./star_bam/`     |
| count      | `./featurecounts/`, `./counts/` |

