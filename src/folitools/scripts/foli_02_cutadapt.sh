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
ulimit -n 1000000

############################################
# User-configurable paths
############################################
FASTP_DIR=./fastp  # input
REST_DIR="./rest_all"
REST_NONDIMER_DIR="./rest"
THREADS=16

GLOB_PATTERN="${1:-*_1.fq.gz}"  # Default if not provided
ADAPTER_DIR="${2:-./data}"

############################################
# Create output directories if needed
############################################
mkdir -p $REST_DIR $REST_NONDIMER_DIR

############################################
# Main loop over paired-end FASTQs in ./fastp
############################################
# UCLA IBD panel
fqr1s=$(ls "$FASTP_DIR"/$GLOB_PATTERN)
i=0
for fqR1 in $fqr1s; do
    i=$((i+1))
    # if [[ $i -le 49 ]]; then
    #     continue
    # fi
    # Derive the matching R2
    fqR2="${fqR1/_1/_2}"
    baseR1=$(basename "$fqR1" .fq.gz)
    sample_name="${baseR1%%_*}"

    # Skip if no matching R2
    if [[ ! -f "$fqR2" ]]; then
        echo "WARNING: Could not find R2 for: $fqR1"
        continue
    fi
    echo "Processing sample: $sample_name"

    cutadapt \
        -j 8 \
        -e 2 \
        -g "file:$ADAPTER_DIR/i5_short.fasta;min_overlap=20" \
        -G "file:$ADAPTER_DIR/i7_short.fasta;min_overlap=20" \
        --rename '{id} {adapter_name} {match_sequence}' \
        --action=none \
        --interleaved \
        -o - \
        "$fqR1" "$fqR2" \
        2> /dev/null \
        | python add_umi.py \
            --o1 "$REST_DIR/${sample_name}_1.fq.gz" \
            --o2 "$REST_DIR/${sample_name}_2.fq.gz"

    # Remove reads that are too short (these are considered primer dimers)
    cutadapt \
    -j 8 \
    --minimum-length 60:60 \
    -o >(paste - - - - | grep -v 'no_adapter' | tr '\t' '\n' | gzip > "$REST_NONDIMER_DIR/${sample_name}_1.fq.gz") \
    -p >(paste - - - - | grep -v 'no_adapter' | tr '\t' '\n' | gzip > "$REST_NONDIMER_DIR/${sample_name}_2.fq.gz") \
    "$REST_DIR/${sample_name}_1.fq.gz" "$REST_DIR/${sample_name}_2.fq.gz" \
    2> /dev/null > /dev/null
done | tqdm --total $(echo "$fqr1s" | wc -w)

if [[ -f "$REST_DIR.stats" ]]; then
    echo "Output file '$REST_DIR.stats' already exists. Skipping seqkit stats."
else
    seqkit stats --all --tabular --threads "$THREADS" $REST_DIR/*.fq.gz > $REST_DIR.stats
fi
if [[ -f "$REST_NONDIMER_DIR.stats" ]]; then
    echo "Output file '$REST_NONDIMER_DIR.stats' already exists. Skipping seqkit stats."
else
    seqkit stats --all --tabular --threads "$THREADS" $REST_NONDIMER_DIR/*.fq.gz > $REST_NONDIMER_DIR.stats
fi

echo "Pre-STAR pipeline completed successfully!"
