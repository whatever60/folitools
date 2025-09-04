#!/usr/bin/env bash

set -euo pipefail

script_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
source "$script_dir/utils.sh"

# Help function
show_help() {
    echo "Usage: foli_04_count.sh <input_files> [output_dir] [threads] [skip]"
    echo "  input_files : Space-separated list of alignment file paths (.bam/.sam)"
    echo "  output_dir  : Output directory for count results (default: ./counts)"
    echo "  threads     : Number of threads to use (default: 16)"
    echo "  skip        : Skip certain steps (default: 0)"
    echo "  -h, --help  : Show this help message"
    exit 0
}

# Check for help option
if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    show_help
fi

INPUT_FILES="${1}"  # Space-separated list of actual file paths
OUTPUT_DIR="${2:-./counts}"  # Output directory for count results
THREADS="${3:-16}"
SKIP="${4:-0}"

COUNTS_DIR="$OUTPUT_DIR"

############################################
# Create output directories if needed
############################################
mkdir -p "$COUNTS_DIR"

# Convert space-separated string back to array
read -ra bams <<< "$INPUT_FILES"

# Validate that all files are alignment format
validate_file_formats "$INPUT_FILES" "alignment"

i=0
for bam in "${bams[@]}"; do
    ((++i))
    if [[ $i -le $SKIP ]]; then
        echo "Skipping sample $i: $bam"
        continue
    fi
    # Extract sample name using utility function
    sample_name="$(extract_sample_name "$bam")"
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
    # umi_tools paired end behaviors:
    #    1. It gets cell tags, gene tags and umi tags only from R1. It does not care if 
    #    R2 has them or not. Also, it won't verify if R1 and R2 have the same tags.
    #    2. Thus, if both R1 and R2 have UMI sequences, users are responsible for concatenating
    #    them into R1 UMI before running umi_tools.
    # umi_tools multimapping behaviors:
    #    1. umi_tools only uses the primary alignment for UMI deduplication and counting.
    #    But secondary/supplementary alignments are also output to the tsv or bam as is.
    #    Therefore, it should be noted that non-primary alignments are not used by umi_tools,
    #    and even they are reflected in the outputs, we do not consider them informative.

    # Run the `group` subcommand to get richer information.
    # Count matrix can be derived from {sample_name}.group.tsv.gz.
    umi_tools group \
        --method unique \
        --per-cell \
        --per-gene \
        --cell-tag CB \
        --cell-tag-split "" \
        --gene-tag XF \
        --umi-tag UC \
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
done | tqdm --total ${#bams[@]} > /dev/null

echo "Pipeline completed successfully!"
