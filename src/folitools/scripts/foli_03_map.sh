REST_DIR="./rest"
STAR_DIR="./star"
BAM_DIR="./star_bam"
THREADS=16
# STAR_INDEX="$HOME/data/gencode/Gencode_human/release_46/STAR_2.7.11b_150"

GLOB_PATTERN="${1:-*_1.fq.gz}"       # Optional pattern for FASTQ files
STAR_INDEX="${2:?Usage: $0 [FASTQ_PATTERN] STAR_INDEX_PATH}"  # Required STAR index path

mkdir -p $STAR_DIR $BAM_DIR

STAR --runThreadN $THREADS \
    --genomeDir $STAR_INDEX \
    --genomeLoad LoadAndExit \
    --outFileNamePrefix _temp/
rm -rf _temp

fqr1s=$(ls "$REST_DIR"/$GLOB_PATTERN)
i=0
for fqR1 in "$REST_DIR"/*_1.fq.gz; do
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


    STAR \
        --runThreadN $THREADS \
        --genomeDir $STAR_INDEX \
        --genomeLoad LoadAndKeep \
        --readFilesIn $fqR1 $fqR2 \
        --readFilesCommand zcat \
        --outFileNamePrefix $STAR_DIR/${sample_name}/ \
        --outFilterMultimapNmax 1000 \
        --outSAMmultNmax 1000 \
        --outSAMunmapped Within \
        --chimOutType WithinBAM \
        --outSAMmode Full \
        --outSAMtype SAM \
        --outSAMorder PairedKeepInputOrder \
        --outStd SAM \
        --outTmpKeep None \
        --quantMode GeneCounts | \
        python add_tags.py \
            --output $BAM_DIR/${sample_name}.bam \
            --cell_tag ${sample_name}
done | tqdm --total $(echo "$fqr1s" | wc -w)

# Unload the STAR genome index to free up memory
STAR \
    --runThreadN "$THREADS" \
    --genomeDir "$STAR_INDEX" \
    --genomeLoad Remove \
    --outFileNamePrefix _temp/
rm -rf _temp

