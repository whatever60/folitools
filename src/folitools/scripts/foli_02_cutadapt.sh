#!/usr/bin/env bash

set -euo pipefail
ulimit -n 1000000

script_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
source "$script_dir/utils.sh"

# Help function
show_help() {
    echo "Usage: foli_02_cutadapt.sh <input_files> [output_dir] <i5_file> <i7_file> [threads] [skip] [delete]"
    echo "  input_files : Space-separated list of FASTQ file paths (.fq/.fastq/.fq.gz/.fastq.gz)"
    echo "  output_dir  : Output directory for UMI-tagged files (default: ./rest_all)"
    echo "  i5_file     : Path to I5 index file"
    echo "  i7_file     : Path to I7 index file"
    echo "  threads     : Number of threads to use (default: 16)"
    echo "  skip        : Skip certain steps (default: 0)"
    echo "  delete      : Delete intermediate files (default: false)"
    echo "  -h, --help  : Show this help message"
    exit 0
}

# Check for help option
if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    show_help
fi

INPUT_FILES="${1}"  # Space-separated list of actual file paths
OUTPUT_DIR="${2:-./rest_all}"  # Output directory for UMI-tagged files
I5_FILE="${3}"
I7_FILE="${4}"
THREADS="${5:-16}"
SKIP="${6:-0}"
DELETE="${7:-false}"

REST_DIR="$OUTPUT_DIR"

mkdir -p $REST_DIR

# Convert space-separated string back to array
read -ra fqr1s <<< "$INPUT_FILES"

# Validate that all files are FASTQ format
validate_file_formats "$INPUT_FILES" "fastq"

i=0
for fqR1 in "${fqr1s[@]}"; do
    ((++i))
    if [[ $i -le $SKIP ]]; then
        echo "Skipping sample $i: $fqR1"
        continue
    fi
    # Derive the matching R2 using utility function
    fqR2="$(derive_r2_from_r1 "$fqR1")"
    sample_name="$(extract_sample_name "$fqR1")"

    # Skip if no matching R2
    if [[ ! -f "$fqR2" ]]; then
        echo "WARNING: Could not find R2 for: $fqR1"
        continue
    fi
    echo "Processing sample: $sample_name"

    cutadapt \
        -j "$THREADS" \
        -e 0.1 \
        -g "file:$I5_FILE;min_overlap=20" \
        -G "file:$I7_FILE;min_overlap=20" \
        --rename '{id} {adapter_name} {match_sequence}' \
        --action none \
        --interleaved \
        -o - \
        "$fqR1" "$fqR2" \
        2> /dev/null \
        | python -m folitools.add_umi \
            --o1 "$REST_DIR/${sample_name}_1.fq.gz" \
            --o2 "$REST_DIR/${sample_name}_2.fq.gz"

    # # Remove reads that are too short or have name containing "no_adapter" (these are 
    # # considered primer dimers)
    # cutadapt \
    #     -j "$THREADS" \
    #     --minimum-length 60:60 \
    #     -o \
    #         >(paste - - - - \
    #         | grep -v 'no_adapter' \
    #         | tr '\t' '\n' \
    #         | gzip \
    #         > "$REST_NONDIMER_DIR/${sample_name}_1.fq.gz") \
    #     -p \
    #         >(paste - - - - \
    #         | grep -v 'no_adapter' \
    #         | tr '\t' '\n' \
    #         | gzip \
    #         > "$REST_NONDIMER_DIR/${sample_name}_2.fq.gz") \
    #     "$REST_DIR/${sample_name}_1.fq.gz" "$REST_DIR/${sample_name}_2.fq.gz" \
    #     &> /dev/null

    # Remove input FASTQs to save space
    if [[ "$DELETE" == "True" ]]; then
        rm "$fqR1" "$fqR2"
    fi

done | tqdm --total ${#fqr1s[@]} > /dev/null

run_seqkit_stats "$REST_DIR"
# run_seqkit_stats "$REST_NONDIMER_DIR"

echo "Pre-STAR pipeline completed successfully!"
