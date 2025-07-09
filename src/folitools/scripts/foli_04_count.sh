#!/usr/bin/env bash

set -euo pipefail

BAM_DIR="./star_bam"  # input
FEATURECOUNTS_DIR="./featurecounts"
COUNTS_DIR="./counts"
# GTF="$HOME/data/gencode/Gencode_mouse/release_M32/gencode.vM32.primary_assembly.annotation.gtf.gz"
# GTF="$HOME/data/gencode/Gencode_human/release_46/gencode.v46.primary_assembly.annotation.gtf.gz"
THREADS=16

GLOB_PATTERN="${1:-*.bam}"                               # Optional BAM glob
GTF="${2:?Usage: $0 [BAM_PATTERN] GTF_PATH}"             # Required GTF file path

############################################
# Create output directories if needed
############################################
mkdir -p "$FEATURECOUNTS_DIR" "$COUNTS_DIR"

bams=$(ls "$BAM_DIR"/$GLOB_PATTERN)
i=0
for star_bam in $bams; do
    ((++i))
    # if [[ $i -le 49 ]]; then
    #     continue
    # fi
    base_name=$(basename $star_bam .bam)
    sample_name="${base_name%%_*}"
    echo "Processing sample: $sample_name"

    # -p: sequencing data is paired-end
    # -B: Only count reads that are properly paired
    # -C: Do not count read pairs that have two ends mapping to different chromosomes or 
    # mapping to the same chromosome but on different strands.
    featureCounts \
        -T "$THREADS" \
        -a "$GTF" \
        -o "$FEATURECOUNTS_DIR/${sample_name}" \
        -p -B -C \
        -R BAM \
        -Rpath "$FEATURECOUNTS_DIR/${sample_name}.bam" \
        "$BAM_DIR/${sample_name}.bam" 2> $FEATURECOUNTS_DIR/$sample_name.log
    # Sort BAM files from featureCounts to prepare for UMI-tools
    sambamba sort "$FEATURECOUNTS_DIR/${sample_name}.bam" 2> /dev/null
    rm "$FEATURECOUNTS_DIR/${sample_name}.bam"

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
        -S "$COUNTS_DIR/${sample_name}.tsv.gz" > $COUNTS_DIR/${sample_name}.log

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
