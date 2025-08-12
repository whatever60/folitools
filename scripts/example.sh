
#!/bin/bash

# Example folitools pipeline with flexible input paths
# Set these environment variables before running
# export STAR_INDEX="/path/to/star/index"
# export GTF_FILE="/path/to/annotation.gtf"

# Step 1: Quality control and preprocessing
# Input: Raw FASTQ files (specify R1 files, R2 files are automatically detected)
foli qc \
    --input "./fastq/*_R1_001.fastq.gz" \
    --output-dir "./fastp" \
    --cores 16

# Step 2: Probe assignment and UMI extraction
# Input: Trimmed FASTQ files from step 1
foli assign-probes \
    --input "./fastp/*_1.fq.gz" \
    --output-dir "./rest_all" \
    --probe-dir ./data/probe \
    --cores 16

# Step 3: Get read statistics (optional)
# Input: UMI-tagged FASTQ files from step 2
foli get-read-stats \
    --input "./rest_all/*_1.fq.gz" \
    --output-dir "./rest_all_stats" \
    --cores 16 \
    --overwrite

# Step 4: Mapping and feature counting
# Input: UMI-tagged FASTQ files from step 2
foli map \
    --input "./rest_all/*_1.fq.gz" \
    --output-dir "./featurecounts" \
    --star-index "$STAR_INDEX" \
    --gtf "$GTF_FILE" \
    --cores 16

# Step 5: UMI-based gene counting
# Input: Sorted BAM files from step 4
foli count \
    --input "./featurecounts/*.sorted.bam" \
    --output-dir "./counts" \
    --cores 16

# Step 6: Generate final count matrix
# Input: UMI grouping files from step 5
foli get-count-mtx \
    --input "./counts/*.group.tsv.gz" \
    --output "./foli_counts.tsv" \
    --gtf "$GTF_FILE"
