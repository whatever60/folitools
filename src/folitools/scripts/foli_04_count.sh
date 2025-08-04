#!/usr/bin/env bash

set -euo pipefail

FEATURECOUNTS_DIR="./featurecounts"  # input
COUNTS_DIR="./counts"

GLOB_PATTERN="${1:-*.bam}"
THREADS="${2:-16}"

############################################
# Create output directories if needed
############################################
mkdir -p "$COUNTS_DIR"
bams=$(ls "$FEATURECOUNTS_DIR"/$GLOB_PATTERN)
i=0
for bam in $bams; do
    ((++i))
    # if [[ $i -le 49 ]]; then
    #     continue
    # fi
    # split by dot to get the sample name
    base_bam=$(basename "$bam" .bam)
    sample_name="${base_bam%%.*}"
    echo "Processing sample: $sample_name"

    # umi_tools behaviors:
    # umi_tools group can take BAM from stdin and output BAM to stdout.
    # umi_tools count can take BAM from stdin and output TSV to stdout, but it does not output BAM.
    # --no-sort-output: umi_tools resort input BAM by read start position (alignment start - soft-clipping) 
    # to as an intermediate step. By default it will sort again by alignment start for output, 
    # and this flag prevents the second sort.
    # --unmapped-reads [discard (default)|use|output]
    # --chimerric-pairs [discard|use (default)|output]
    # --unpaired-reads [discard|use (default)|output]
    # --log <file> will append to the file instead of overwriting it. So we use 
    # --log2stderr with 2> <file> to redirect stderr to a clean file.
    # What is the effect of --paired?
    # What is the effect of --out-sam in umi_tools count if it does not output BAM at all?

    # Run the `group` subcommand to get richer information.
    # Count matrix can be derived from {sample_name}.group.tsv.gz.
    umi_tools group \
        --method unique \
        --per-cell \
        --per-gene \
        --cell-tag CB \
        --cell-tag-split "" \
        --gene-tag XT \
        --umi-tag UC \
        --assigned-status-tag XT \
        --extract-umi-method tag \
        --paired \
        --unmapped-reads output \
        --chimeric-pairs output \
        --unpaired-reads output \
        --no-sort-output \
        --group-out "$COUNTS_DIR/${sample_name}.group.tsv.gz" \
        --stdin $bam \
        --log2stderr \
        2> "$COUNTS_DIR/${sample_name}.group.log"
        # --skip-tags-regex "(?!)" \
        # | sambamba sort \
        #     --nthreads "$((THREADS-1))" \
        #     --memory-limit 16GB \
        #     --out "$COUNTS_DIR/${sample_name}.sorted.bam" \
        #     /dev/stdin
done | tqdm --total $(echo "$bams" | wc -w) > /dev/null

echo "Pipeline completed successfully!"
