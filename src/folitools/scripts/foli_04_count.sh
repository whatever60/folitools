#!/usr/bin/env bash

set -euo pipefail

BAM_DIR="./star_bam"  # input
FEATURECOUNTS_DIR="./featurecounts"
COUNTS_DIR="./counts"

GLOB_PATTERN="${1:-*.bam}"
THREADS="${2:-16}"

############################################
# Create output directories if needed
############################################
mkdir -p "$COUNTS_DIR"

bams=$(ls "$FEATURECOUNTS_DIR"/$GLOB_PATTERN)
i=0
for star_bam in $bams; do
    ((++i))
    # if [[ $i -le 49 ]]; then
    #     continue
    # fi
    # split by dot to get the sample name
    base_bam=$(basename "$star_bam" .bam)
    sample_name="${base_bam%%.*}"
    echo "Processing sample: $sample_name"

    umi_tools count \
        --method unique \
        --per-cell \
        --per-gene \
        --cell-tag CB \
        --cell-tag-split "" \
        --gene-tag XT \
        --umi-tag UP \
        --assigned-status-tag XS \
        --wide-format-cell-counts \
        --extract-umi-method tag \
        -I "$FEATURECOUNTS_DIR/${sample_name}.sorted.bam" \
        -S "$COUNTS_DIR/${sample_name}.tsv.gz" \
        > $COUNTS_DIR/${sample_name}.log

    # Run the `group` subcommand to get richer information.
    # Output: grouped.bam and group.tsv.gz
    # Count matrix can be derived from group.tsv.gz.
    # So the above `count` command becomes optional if you run `group`. But if you don't need 
    # per-read information, you just need to run `count`.
    umi_tools group \
        --method=unique \
        --per-cell \
        --per-gene \
        --cell-tag CB \
        --cell-tag-split "" \
        --gene-tag XT \
        --umi-tag UP \
        --assigned-status-tag XS \
        --extract-umi-method tag \
        --group-out "$COUNTS_DIR/${sample_name}.group.tsv.gz" \
        --output-bam \
        -S "$COUNTS_DIR/${sample_name}.grouped.bam" \
        -I "$FEATURECOUNTS_DIR/${sample_name}.sorted.bam" \
        > "$COUNTS_DIR/${sample_name}.group.log"

done | tqdm --total $(echo "$bams" | wc -w) > /dev/null

echo "Pipeline completed successfully!"
