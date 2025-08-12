#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$script_dir/utils.sh"

INPUT_FILES="${1}"  # Space-separated list of actual file paths
OUTPUT_DIR="${2:-./fastp}"  # Output directory for trimmed files
THREADS="${3:-16}"

FASTP_DIR="$OUTPUT_DIR"
FASTP_FASTQC_DIR="${OUTPUT_DIR}_fastqc"

mkdir -p "$FASTP_DIR" "$FASTP_FASTQC_DIR"

# Create a temporary directory for input FastQC (only if files are found)
FASTQ_FASTQC_DIR="./fastq_fastqc"
mkdir -p "$FASTQ_FASTQC_DIR"

# Convert space-separated string back to array
read -ra fqr1s <<< "$INPUT_FILES"

# Generate stats for input files from their directories (run once per unique directory)
input_dirs=$(printf '%s\n' "${fqr1s[@]}" | xargs -I {} dirname {} | sort -u)
for dir in $input_dirs; do
    if [[ -d "$dir" ]] && ls "$dir"/*.{fastq.gz,fq.gz,fastq,fq} 1> /dev/null 2>&1; then
        run_seqkit_stats "$dir"
    fi
done

############################################
# Main loop over paired-end FASTQs
############################################
i=0
for fqR1 in "${fqr1s[@]}"; do
    ((++i))
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
    # if [[ -f "$trimmed_R1" && -f "$trimmed_R2" ]]; then
    #     echo "Skipping already processed sample: $sample_name"
    #     continue
    # fi

    # Skip if no matching R2
    if [[ ! -f "$fqR2" ]]; then
        echo "WARNING: Could not find R2 for: $fqR1"
        continue
    fi
    echo "Processing sample: $sample_name"

    # run_fastqc "$fqR1" "$fqR2" "$FASTP_FASTQC_DIR" "$THREADS"

    # Run fastqc "$fqR1" "$fqR2" "$FASTQ_FASTQC_DIR" "$THREADS"
    # No need for poly-G or TruSeq adapter trimming, fastp trimming by overlapp analysis is good enough. 
    # fastp \
    #     --in1 "$fqR1" \
    #     --in2 "$fqR2" \
    #     --stdout \
    #     --thread "$THREADS" \
    #     --cut_tail \
    #     --correction \
    #     --html "$FASTP_DIR/${sample_name}.fastp.html" \
    #     --json "$FASTP_DIR/${sample_name}.fastp.json" \
    #     2> /dev/null | cutadapt \
    #         -e 2 \
    #         -a "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=5" \
    #         -A "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT;min_overlap=5" \
    #         -j "$THREADS" \
    #         --interleaved \
    #         -o "$trimmed_R1" \
    #         -p "$trimmed_R2" \
    #         - &> /dev/null
        # --trim_poly_g \
        # --trim_poly_x \
        # --adapter_sequence AGATCGGAAGAGCACACGTC \
        # --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGT \
        # Adapter sequences are the first 20bp of Illumina TruSeq adapters:
        # i7 adapter: AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
        # i5 adapter: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

    fastp \
        --in1 "$fqR1" \
        --in2 "$fqR2" \
        --out1 "$trimmed_R1" \
        --out2 "$trimmed_R2" \
        --thread "$THREADS" \
        --cut_tail \
        --correction \
        --html "$FASTP_DIR/${sample_name}.fastp.html" \
        --json "$FASTP_DIR/${sample_name}.fastp.json" \
        2> /dev/null

    run_fastqc "$trimmed_R1" "$trimmed_R2" "$FASTP_FASTQC_DIR" "$THREADS"

    # Remove input FASTQs to save space
    # rm "$fqR1" "$fqR2"
done | tqdm --total ${#fqr1s[@]} > /dev/null

run_seqkit_stats "$FASTP_DIR"

echo "Preprocessing pipeline completed successfully!"
