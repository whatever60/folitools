#!/usr/bin/env bash

###############################################################################
# RNA-seq Pre-STAR Pipeline
#
# Steps:
# 1. FastQC + SeqKit on raw reads
# 2. fastp trimming (3' sliding window)
# 3. FastQC + SeqKit on trimmed reads
# 4. Cutadapt for demultiplexing and UMI extraction
# 5. Combine reads for STAR alignment
#
###############################################################################

# set -euo pipefail

############################################
# User-configurable paths
############################################
FASTQ_DIR="./fastq"
FASTQ_FASTQC_DIR="./fastq_fastqc"
FASTP_DIR="./fastp"
FASTP_FASTQC_DIR="./fastp_fastqc"
THREADS=8

GLOB_PATTERN="${1:-*_R1_*.fastq.gz}"  # Use first argument if given, otherwise default

############################################
# Create output directories if needed
############################################
mkdir -p "$FASTQ_DIR" "$FASTQ_FASTQC_DIR" "$FASTP_DIR" "$FASTP_FASTQC_DIR"

############################################
# SeqKit stats for raw FASTQs
############################################
if [[ -f "$FASTQ_DIR.stats" ]]; then
    echo "Output file '$FASTQ_DIR.stats' already exists. Skipping seqkit stats."
else
    seqkit stats --all --tabular --threads "$THREADS" "$FASTQ_DIR"/*.fastq.gz > $FASTQ_DIR.stats
fi

############################################
# Main loop over paired-end FASTQs in ./fastq
############################################
fqr1s=$(ls "$FASTQ_DIR"/$GLOB_PATTERN)
i=0
for fqR1 in $fqr1s; do
    i=$((i+1))
    # if [[ $i -le 49 ]]; then
    #     continue
    # fi
    
    # Derive the matching R2
    fqR2="${fqR1/_R1_/_R2_}"
    baseR1=$(basename "$fqR1" .fastq.gz)
    sample_name="${baseR1%%_*}"
    trimmed_R1="$FASTP_DIR/${sample_name}_1.fq.gz"
    trimmed_R2="$FASTP_DIR/${sample_name}_2.fq.gz"

    # skip if "$FASTP_DIR/${sample_name}_1.fq.gz" already exists
    if [[ -f "$trimmed_R1" && -f "$trimmed_R2" ]]; then
        echo "Skipping already processed sample: $sample_name"
        continue
    fi

    # Skip if no matching R2
    if [[ ! -f "$fqR2" ]]; then
        echo "WARNING: Could not find R2 for: $fqR1"
        continue
    fi
    echo "Processing sample: $sample_name"

    # ########################################################
    # # Step 1: FastQC on raw reads
    # ########################################################
    seqkit_out=$(seqkit seq $fqR1 | head -n 1 | wc -l)
    if [[ $seqkit_out -gt 0 ]]; then
        fastqc -t "$THREADS" -o "$FASTQ_FASTQC_DIR" "$fqR1" "$fqR2" 2> /dev/null
    fi
    ########################################################
    # Step 2: fastp for trimming
    ########################################################
    # We use novaseq, so we trim poly G tail
    fastp \
        --in1 "$fqR1" \
        --in2 "$fqR2" \
        --stdout \
        --thread "$THREADS" \
        --cut_tail \
        --correction \
        --html "$FASTP_DIR/${sample_name}.fastp.html" \
        --json "$FASTP_DIR/${sample_name}.fastp.json" \
        2> /dev/null | cutadapt \
            -e 2 \
            -a "AGATCGGAAGAGCACACGTC;min_overlap=5" \
            -A "AGATCGGAAGAGCGTCGTGT;min_overlap=5" \
            -j "$THREADS" \
            --interleaved \
            -o "$trimmed_R1" \
            -p "$trimmed_R2" \
            - &> /dev/null
        # --trim_poly_g \
        # --trim_poly_x \
        # --adapter_sequence AGATCGGAAGAGCACACGTC \
        # --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGT \

        # Adapter sequences are the first 20bp of Illumina TruSeq adapters:
        # i7 adapter: AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
        # i5 adapter: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

    # ########################################################
    # # Step 3: FastQC on trimmed reads
    # ########################################################
    seqkit_out=$(seqkit seq $trimmed_R1 | head -n 1 | wc -l)
    if [[ $seqkit_out -gt 0 ]]; then
        fastqc -t "$THREADS" -o "$FASTP_FASTQC_DIR" "$trimmed_R1" "$trimmed_R2" 2> /dev/null
    fi
done | tqdm --total $(echo "$fqr1s" | wc -w)

if [[ -f "$FASTP_DIR.stats" ]]; then
    echo "Output file '$FASTP_DIR.stats' already exists. Skipping seqkit stats."
else
    seqkit stats --all --tabular --threads "$THREADS" "$FASTP_DIR"/*.fq.gz > $FASTP_DIR.stats
fi

echo "Pre-STAR pipeline completed successfully!"
