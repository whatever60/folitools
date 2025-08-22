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
foli qc \
    --input "/path/to/fastq/*_R1_*.fastq.gz" \
    --output-dir "/path/to/trimmed_reads"
```

This command performs read trimming using `fastp` and runs quality check on output files with `fastqc` and `seqkit`.

You can specify input FASTQ files from any location using absolute paths, relative paths, or glob patterns. The paired-end R2 files are automatically derived from R1 files by replacing `_R1_` with `_R2_`.

Output files:
- Trimmed reads in the specified output directory
- FastQC reports for trimmed reads
- Summary statistics for trimmed reads

### Step 2. Probe assignment

```bash
foli assign-probes \
    --input "/path/to/trimmed_reads/*_1.fq.gz" \
    --output-dir "/path/to/umi_tagged" \
    --i5 /path/to/i5_short.fasta \
    --i7 /path/to/i7_short.fasta
```

This command performs gene probe primer assignment and trimming using `cutadapt`. This step extracts UMI sequences from adapter matches and embeds them in read names for downstream processing. Quality statistics are also generated for the output files using `seqkit`.

After probe assignment, you can get read statistics:

```bash
foli get-read-stats \
    --input "/path/to/umi_tagged/*_1.fq.gz" \
    --output-dir "/path/to/read_stats"
```

This analyzes the UMI-tagged reads to generate detailed statistics about primer assignment and read characteristics.

Output files:
- UMI-tagged, adapter-trimmed reads
- Summary statistics for UMI-tagged reads
- Read statistics in parquet format

### Step 3. Mapping and Feature Counting

```bash
foli map \
    --input "/path/to/umi_tagged/*_1.fq.gz" \
    --output-bam "/path/to/mapped_reads" \
    --output-star "/path/to/star" \
    --star-index /path/to/star_index \
    --gtf /path/to/annotation.gtf
```

This command performs streamlined alignment and feature counting using `STAR` and `featureCounts`. This step filters out short reads (< 60bp, considered primer dimers) using `cutadapt`, aligns filtered reads to the genome using `STAR`, and assigns reads to genomic features using `featureCounts`.

For fractional counting of reads that overlap multiple features, you can use the `--allow-overlap` and `--allow-multimapping` flags:
```bash
foli map \
    --input "/path/to/umi_tagged/*_1.fq.gz" \
    --output-bam "/path/to/mapped_reads" \
    --output-star "/path/to/star" \
    --star-index /path/to/star_index \
    --gtf /path/to/annotation.gtf \
    --allow-overlap \
    --allow-multimapping
```

UMI sequences and cell barcodes are added to BAM records as tags (`US` for raw UMI, `UC` for filtered UMI, `CB` for cell barcode) for the next counting step.

Output files:
- STAR alignment files in the specified STAR output directory
- Sorted, tagged BAM files in the specified BAM output directory
- `featureCounts` tables and logs

### Step 4. UMI-based Gene Counting

```bash
foli count \
    --input "/path/to/mapped_reads/*.sorted.bam" \
    --output-dir "/path/to/count_results"
```

This command processes BAM files using `umi_tools group` to generate UMI-deduplicated count data. This step:

1. Groups reads by UMI, cell barcode, and gene assignment
2. Handles paired-end reads, unmapped reads, and chimeric pairs
3. Outputs detailed grouping information for downstream analysis

Input BAM files should contain UMI tags (`UC`), cell tags (`CB`), and gene assignment tags (`XT`) from the previous step.

After counting, generate the final count matrix:

```bash
foli get-count-mtx \
    --input "/path/to/count_results/*.group.tsv.gz" \
    --output "/path/to/final_counts.tsv" \
    --gtf /path/to/annotation.gtf
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

## Primer Selection Functionality

The primer selection module provides tools for designing and recovering PCR primer sets for targeted amplicon sequencing experiments.

### Workflow

The primer selection workflow provides a complete pipeline for primer design:

```bash
# Run complete primer design workflow using built-in reference
foli-primer workflow \
    --genes genes.tsv \
    --species mouse \
    --output-dir primer_design/

# Or use custom transcriptome FASTA for better coverage
foli-primer workflow \
    --genes genes.tsv \
    --txome-fasta /path/to/transcriptome.fasta \
    --output-dir primer_design/
```

**Input**: Gene table TSV with columns `gene` and `group`  
**Output**: Complete primer design including optimized primer sets, amplicon sequences, and IDT-compatible ordering files

**Note**: You can use either `--species` (mouse/human) for built-in references or `--txome-fasta` for custom transcriptome files. Custom FASTA files often provide better gene coverage than the built-in Gencode references.

### Recover

The recover functionality helps validate and analyze primer sets from IDT order files:

```bash
# Recover primer information using built-in reference
foli-primer recover \
    --order-excel idt_order.xlsx \
    --output-dir recovered_output/ \
    --species human

# Or use custom transcriptome FASTA (recommended for better coverage)
foli-primer recover \
    --order-excel idt_order.xlsx \
    --output-dir recovered_output/ \
    --txome-fasta /path/to/transcriptome.fasta
```

This command analyzes primer sequences from an IDT order file, validates them against a reference transcriptome, and generates output files for downstream analysis.

**Note**: Using `--txome-fasta` with a custom transcriptome file is recommended over the built-in `--species` references as it typically provides better gene coverage than the packaged Gencode references.

You can also specify the amplicon length range, and whether the primers have linkers:
```bash
foli-primer recover \
    --order-excel idt_order.xlsx \
    --output-dir recovered_output/ \
    --txome-fasta /path/to/transcriptome.fasta \
    --has-linker \
    --amplicon-length-range 300 400
```

**Input**: IDT order Excel file with primer sequences  
**Output**: 
- `summary_primer_to_order.xlsx`: Validated primer summary with amplicon information
- `primer_diagnose.pdf`: PDF report with analysis plots and statistics
- `i5_short.fasta` and `i7_short.fasta`: FASTA files for sequencing adapters

### Dimer Evaluation (through Python)

A primer set can be evaluated by calculating the thermodynamic properties of potential primer dimers for each pair of primers.

```python
from folitools.primer_selection.eval_dimer import dimer_thermo_property

# Evaluate primer-dimer interactions
result_matrix = dimer_thermo_property(
    primer_fwd_fasta="forward_primers.fasta",
    primer_rev_fasta="reverse_primers.fasta", 
    output_dir="dimer_analysis/",
    output_suffix="_analysis"
)
```

**Input**: FASTA files containing forward and reverse primer sequences  
**Output**: 
- Pairwise interaction matrix with thermodynamic properties
- Detailed analysis files in the output directory

## Testing

Run all tests:
```bash
pytest
```

Run only shell script tests:
```bash
pytest tests/test_shell_scripts.py -v
```

Run tests with coverage:
```bash
pytest --cov=src/folitools --cov-report=html
```

## TODO

- Add QC report code that summarizes all metrics into a table (need refactoring)
