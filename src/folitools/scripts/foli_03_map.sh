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
REST_DIR="./rest_all"
STAR_DIR="./star"
FEATURECOUNTS_DIR="./featurecounts"

echo "Creating output directories..."
mkdir -p "$STAR_DIR" "$FEATURECOUNTS_DIR"

# --- Dynamic Core Allocation ---
# read STAR_THREADS FC_THREADS < <(python -m folitools.scripts.foli_03_map_utils --total-cores "$((CORES - 2))")
# echo "Core allocation -> STAR: $STAR_THREADS; featureCounts: $FC_THREADS"

# --- Pre-load STAR Genome Index ---
echo "Loading STAR genome into memory..."
STAR --runThreadN "$STAR_THREADS" \
    --genomeDir "$STAR_INDEX" \
    --genomeLoad LoadAndExit \
    --outFileNamePrefix _temp/ \
    > /dev/null
rm -rf _temp

fqr1s=$(ls "$REST_DIR"/$GLOB_PATTERN)
i=0
for fqR1 in $fqr1s; do
    ((++i))
    # if [[ $i -le 7 ]]; then
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

    # Lessons learned during crafting the code below:
    # The big command below basically drops short reads using cutadapt and directly pipes its 
    # output to STAR, without saving any intermediate files.
    # FIFO files are used for communication between processes. To make it work, I had to
    # 1. Explicitly say `exec 3<>"$FIFO_R1"` and `exec 4<>"$FIFO_R2"` after creating the FIFOs. This it 
    #  to ensure that the FIFOs are opened for both reading and writing and to prevent blocking. As a
    #  result, we also need to explicitly close the FIFOs after the cutadapt command (normally the FIFO 
    #  will close when it sees EOF), so that STAR knows that the input is complete. Suprisingly, the way 
    #  to do this is to echo a newline to the FIFOs instead of the moore common `exec 3>&- 4>&-`.
    #  I really don't know why.
    # 2. Run the consumer processes (STAR) first in the background (&) and then the producer (cutadapt).
    #   Otherwise, blocking still occurs, though I don't understand why.
    # 3. Record the PID of the STAR process and wait for it to finish before removing the FIFOs and running #   next steps. Otherwise the output of STAR might be incomplete.
    # 4. Somehow cutadapt does not work well with FIFO files. That is, if I set -o and -p to FIFOs, it
    #  will just block. Instead, I let cutadapt output interleaved FASTQ to stdout and used another awk
    #  command to split it back into two FIFOs.
    # 5. FIFO file names must be explicit since STAR reads data based on file name extension. We not only
    #  need the FIFOs to be text (not gzipped), but also need to set --readFilesCommand to cat. Just
    #  setting --readFilesCommand to cat still results in error when the FIFOs are named with .gz.

    FIFO_R1="$FEATURECOUNTS_DIR/${sample_name}_1.fifo.fq"
    FIFO_R2="$FEATURECOUNTS_DIR/${sample_name}_2.fifo.fq"
    # Delete if they already exist
    rm -f "$FIFO_R1" "$FIFO_R2"
    mkfifo "$FIFO_R1" "$FIFO_R2"
    exec 3<>"$FIFO_R1"    # O_RDWR on FIFO_R1
    exec 4<>"$FIFO_R2"    # O_RDWR on FIFO_R2
    # Remove reads that are too short (these are considered primer dimers)
    # About STAR arguments:
    # --outSAMtype BAM or other BAM types are only compatible with --outSAMorder Paired.
    #   So we will leave it as this default.
    # Note that we are sorting afterwards with sambamba. So no need to sort in STAR or featureCounts.
    STAR \
        --runThreadN "$((CORES - 1))" \
        --genomeDir "$STAR_INDEX" \
        --genomeLoad LoadAndKeep \
        --readFilesIn "$FIFO_R1" "$FIFO_R2" \
        --readFilesCommand cat \
        --outFileNamePrefix "$STAR_DIR/${sample_name}/" \
        --outFilterMultimapNmax 1000 \
        --outSAMmultNmax 1000 \
        --outSAMunmapped Within \
        --chimOutType WithinBAM \
        --outSAMmode Full \
        --outSAMtype BAM Unsorted \
        --outSAMorder Paired \
        --outStd BAM_Unsorted \
        --outBAMcompression 1 \
        --outTmpKeep None \
        --quantMode GeneCounts | \
        python -m folitools.add_tags \
            --cell_tag ${sample_name} \
            --output "$FEATURECOUNTS_DIR/${sample_name}.temp.bam" &
    star_pid=$!
    cutadapt \
        -j 1 \
        --interleaved \
        --minimum-length 60:60 \
        "$REST_DIR/${sample_name}_1.fq.gz" "$REST_DIR/${sample_name}_2.fq.gz" \
        2> /dev/null \
        | awk -v p1="$FIFO_R1" -v p2="$FIFO_R2" \
            '
                {
                    if (((NR - 1) % 8) < 4) {
                        print $0 > p1
                    } else {
                        print $0 > p2
                    }
                }
            '
    # Write EOF to the FIFOs to signal completion
    echo -e "\n" >&3
    echo -e "\n" >&4
    # exec 3>&- 4>&-  # Close the FIFOs
    wait "$star_pid"
    rm -f "$FIFO_R1" "$FIFO_R2"

    # About featureCounts arguments:
    # -T: Number of threads
    # -p: sequencing data is paired-end
    # -B: Only count reads that are properly paired
    # -C: Do not count read pairs that have two ends mapping to different chromosomes or 
    #   mapping to the same chromosome but on different strands.
    # Arguments that might be relevant but we leave as default:
    # -s 0: Strand specificity (0 = unstranded, 1 = stranded, 2 = reversely stranded)
    # -t exon: Use exon feature type for counting
    # -g gene_id: Use gene_id attribute for counting
    # Note that setting -R/--Rpath is not allowed when featureCounts reads BAM/SAM from stdin.
    #  Also the BAM sorting step next cannot be streamlined, as --Rpath does not support 
    #  stdout. That's why we have to have some temp BAM files here.
    # featureCounts do not expect the BAM to be ordered in any ways, but it's best for it to have
    #  paired reads close together (featureCounts will try look for read pairs, but not 
    #  exhaustively), which is apprently true for name sorted bam, and also mostly true 
    #  for coordinate sorted BAMs. In the above STAR output, read pairs are always 
    #  together (R2 follows R1), but there is no ordering across different pairs, which does not matter.
    #  By default, featureCounts will sort such that paired reads are together. We here 
    #  turn on the --donotsort flag, as paired reads are already together. Even though it won't make 
    #  any difference here, it's still good to be explicit.
    featureCounts \
        -T "$CORES" \
        -a "$GTF_PATH" \
        -o "$FEATURECOUNTS_DIR/${sample_name}.txt" \
        -p -B -C \
        --donotsort \
        -R BAM \
        --Rpath "$FEATURECOUNTS_DIR/" \
        "$FEATURECOUNTS_DIR/${sample_name}.temp.bam" \
        2> $FEATURECOUNTS_DIR/$sample_name.log

    mv "$FEATURECOUNTS_DIR/${sample_name}.temp.bam.featureCounts.bam" "$FEATURECOUNTS_DIR/${sample_name}.bam"
    sambamba sort \
        --nthreads $CORES \
        --memory-limit 16GB \
        "$FEATURECOUNTS_DIR/${sample_name}.bam" \
        2> /dev/null

    rm "$FEATURECOUNTS_DIR/${sample_name}.temp.bam" "$FEATURECOUNTS_DIR/${sample_name}.bam"
# done
done | tqdm --total $(echo "$fqr1s" | wc -w) > /dev/null

# Unload the STAR genome index to free up memory
echo "Unloading STAR genome from memory."
STAR \
    --runThreadN "$STAR_THREADS" \
    --genomeDir "$STAR_INDEX" \
    --genomeLoad Remove \
    --outFileNamePrefix _temp/ \
    > /dev/null
rm -rf _temp

echo "Pipeline finished successfully."