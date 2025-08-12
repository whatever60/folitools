
#!/bin/bash

# Example folitools pipeline with flexible input paths
# Set these environment variables before running
# export STAR_INDEX="/path/to/star/index"
# export GTF_FILE="/path/to/annotation.gtf"

# Step 0 (optional): Input file QC
mapfile -d '' files < <(printf '%s\0' ./fastq/*_R{1,2}_001.fastq.gz | sort -zV)
seqkit stats --all --tabular --threads 16 "${files[@]}" > ./fastq.stats
fastqc -t 16 -o ./fastq_fastqc "${files[@]}" &> /dev/null

# Step 1: Quality control and preprocessing
# Input: Raw FASTQ files (specify R1 files, R2 files are automatically detected)
foli qc \
    --input "./fastq/*_R1_001.fastq.gz" \
    --output-dir "./fastp" \
    --cores 16 \
    --skip 0

# Step 2: Probe assignment and UMI extraction
# Input: Trimmed FASTQ files from step 1
foli assign-probes \
    --input "./fastp/*_1.fq.gz" \
    --output-dir "./rest_all" \
    --i5 ./data/probe/i5_short.fasta \
    --i7 ./data/probe/i7_short.fasta \
    --cores 16 \
    --skip 0

# Get read statistics for step 2
foli get-read-stats \
    --input "./rest_all/*_1.fq.gz" \
    --output-dir "./rest_all_stats" \
    --cores 16 \
    --overwrite \
    --skip 0

# Step 3: Mapping and feature counting
# Input: UMI-tagged FASTQ files from step 2
foli map \
    --input "./rest_all/*_1.fq.gz" \
    --output-bam "./featurecounts" \
    --output-star "./star" \
    --star-index "$STAR_INDEX" \
    --gtf "$GTF_FILE" \
    --cores 16 \
    --skip 0

# Step 4: UMI-based gene counting
# Input: Sorted BAM files from step 3
foli count \
    --input "./featurecounts/*.sorted.bam" \
    --output-dir "./counts" \
    --cores 16 \
    --skip 0

# Generate final count matrix for step 4
foli get-count-mtx \
    --input "./counts/*.group.tsv.gz" \
    --output "./foli_counts.tsv" \
    --gtf "$GTF_FILE"
