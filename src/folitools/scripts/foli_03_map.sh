#!/usr/bin/env bash
set -euo pipefail

# --- Argument Parsing ---
# A more robust argument parsing loop
CORES=""
STAR_INDEX=""
GTF_PATH=""
GLOB_PATTERN="*_1.fq.gz" # Default pattern

# A simple help message
usage() {
    echo "Usage: $0 --cores <int> --star-index <path> --gtf <path> [--pattern <glob>]"
    echo ""
    echo "  --cores        : Total number of cores to allocate for the pipeline."
    echo "  --star-index   : Path to the STAR genome index directory."
    echo "  --gtf          : Path to the GTF annotation file for featureCounts."
    echo "  --pattern      : Optional glob pattern for R1 FASTQ files (default: '*_1.fq.gz')."
    exit 1
}

# Parse named arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --cores) CORES="$2"; shift ;;
        --star-index) STAR_INDEX="$2"; shift ;;
        --gtf) GTF_PATH="$2"; shift ;;
        --pattern) GLOB_PATTERN="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; usage ;;
    esac
    shift
done

# Check for mandatory arguments
if [[ -z "$CORES" || -z "$STAR_INDEX" || -z "$GTF_PATH" ]]; then
    echo "Error: Missing mandatory arguments."
    usage
fi

# STAR_INDEX="$HOME/data/gencode/Gencode_human/release_46/STAR_2.7.11b_150"
# STAR_INDEX="$HOME/data/gencode/Gencode_mouse/release_M32/STAR_2.7.10b_150"
# GTF="$HOME/data/gencode/Gencode_mouse/release_M32/gencode.vM32.primary_assembly.annotation.gtf.gz"
# GTF="$HOME/data/gencode/Gencode_human/release_46/gencode.v46.primary_assembly.annotation.gtf.gz"

# --- Directory Setup ---
REST_DIR="./rest"
STAR_DIR="./star_bam"
FEATURECOUNTS_DIR="./featurecounts"
COUNTS_DIR="./counts"

echo "Creating output directories..."
mkdir -p "$STAR_DIR" "$FEATURECOUNTS_DIR" "$COUNTS_DIR"

# --- Dynamic Core Allocation ---
read STAR_THREADS FC_THREADS < <(python -m folitools.scripts.foli_03_map_utils --total-cores "$((CORES - 1))")
echo " -> STAR Threads: $STAR_THREADS"
echo " -> featureCounts Threads: $FC_THREADS"

# --- Pre-load STAR Genome Index ---
echo "Loading STAR genome into memory..."
STAR --runThreadN "$STAR_THREADS" \
    --genomeDir "$STAR_INDEX" \
    --genomeLoad LoadAndExit \
    --outFileNamePrefix _temp/
rm -rf _temp

fqr1s=$(ls "$REST_DIR"/$GLOB_PATTERN)
i=0
for fqR1 in $fqr1s; do
    ((++i))
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


    # About featureCounts arguments:
    # -T: Number of threads
    # -p: sequencing data is paired-end
    # -B: Only count reads that are properly paired
    # -C: Do not count read pairs that have two ends mapping to different chromosomes or 
    # mapping to the same chromosome but on different strands.
    # arguments that might be relevant but we leave as default:
    # -s 0: Strand specificity (0 = unstranded, 1 = stranded, 2 = reversely stranded)
    # -t exon: Use exon feature type for counting
    # -g gene_id: Use gene_id attribute for counting
    # Note that we are sorting afterwards with sambamba. So no need to sort in STAR or featureCounts.
    STAR \
        --runThreadN "$STAR_THREADS" \
        --genomeDir "$STAR_INDEX" \
        --genomeLoad LoadAndKeep \
        --readFilesIn "$fqR1" "$fqR2" \
        --readFilesCommand zcat \
        --outFileNamePrefix "$STAR_DIR/${sample_name}/" \
        --outFilterMultimapNmax 1000 \
        --outSAMmultNmax 1000 \
        --outSAMunmapped Within \
        --chimOutType WithinBAM \
        --outSAMmode Full \
        --outSAMtype BAM Unsorted \
        --outSAMorder PairedKeepInputOrder \
        --outBAMcompression 1 \
        --outStd BAM_Unsorted \
        --outTmpKeep None \
        --quantMode GeneCounts | \
        python -m folitools.add_tags \
            --cell_tag ${sample_name} | \
        featureCounts \
            -T "$FC_THREADS" \
            -a "$GTF_PATH" \
            -o "$COUNTS_DIR/${sample_name}.txt" \
            -p -B -C \
            --donotsort \
            -R BAM \
            -Rpath "$FEATURECOUNTS_DIR/${sample_name}.bam" \
            - \
            2> $FEATURECOUNTS_DIR/$sample_name.log

    sambamba sort \
        --nthreads $CORES \
        --memory-limit 16GB \
        "$FEATURECOUNTS_DIR/${sample_name}.bam" \
        2> /dev/null

    rm "$FEATURECOUNTS_DIR/${sample_name}.bam"
done | tqdm --total $(echo "$fqr1s" | wc -w) > /dev/null

# Unload the STAR genome index to free up memory
echo "Unloading STAR genome from memory."
STAR \
    --runThreadN "$STAR_THREADS" \
    --genomeDir "$STAR_INDEX" \
    --genomeLoad Remove \
    --outFileNamePrefix _temp/
rm -rf _temp

echo "Pipeline finished successfully."