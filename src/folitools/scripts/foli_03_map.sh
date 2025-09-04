#!/usr/bin/env bash
set -euo pipefail

script_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
source "$script_dir/utils.sh"

# Help function
show_help() {
    echo "Usage: $0 <input_files> <output_bam_dir> <output_star_dir> <star_index> <gtf_file> [threads] [skip] [delete] [strand] [allow_overlap] [allow_multimapping]"
    echo "  input_files    : Space-separated list of R1 FASTQ (.fq/.fastq/.fq.gz/.fastq.gz) or BAM/SAM file paths"
    echo "  output_bam_dir : Output directory for BAM files"
    echo "  output_star_dir: Output directory for STAR files (required if any FASTQ inputs)"
    echo "  star_index     : Path to the STAR genome index directory (required if any FASTQ inputs)"
    echo "  gtf_file       : Path to the GTF annotation file"
    echo "  threads        : Number of threads (default: 1)"
    echo "  skip           : Number of samples to skip (default: 0)"
    echo "  delete         : Delete input files after processing (default: false)"
    echo "  strand         : Strand specificity for featureCounts (0=unstranded, 1=stranded, 2=reversely stranded) (default: 0)"
    echo "  allow_overlap  : Allow reads to be assigned to overlapping features (default: false)"
    echo "  allow_multimapping : Allow multi-mapping reads to be counted (default: false)"
    echo "  -h, --help     : Show this help message"
    exit 0
}

# Check for help option
if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    show_help
fi

# Parse positional arguments
INPUT_FILES="${1:-}"
OUTPUT_BAM="${2:-}"
OUTPUT_STAR="${3:-}"
STAR_INDEX="${4:-}"
GTF_PATH="${5:-}"
THREADS="${6:-1}"
SKIP="${7:-0}"
DELETE="${8:-false}"
STRAND="${9:-0}"
ALLOW_OVERLAP="${10:-false}"
ALLOW_MULTIMAPPING="${11:-false}"

if [[ -z "$INPUT_FILES" || -z "$OUTPUT_BAM" || -z "$GTF_PATH" ]]; then
    show_help
fi

# Validate strand parameter
if [[ ! "$STRAND" =~ ^[0-2]$ ]]; then
    echo "ERROR: strand parameter must be 0, 1, or 2"
    exit 1
fi

# Convert space-separated string back to array to analyze input files
read -ra input_files <<< "$INPUT_FILES"

# Validate file formats
validate_file_formats "$INPUT_FILES" "both"

# Analyze input files to determine what we need
PRELOAD=false
HAS_FASTQ=false
HAS_BAM=false

for file in "${input_files[@]}"; do
    if is_fastq_file "$file"; then
        HAS_FASTQ=true
    elif is_alignment_file "$file"; then
        HAS_BAM=true
    fi
done

# Determine if we need to preload STAR
if [[ "$HAS_FASTQ" == true ]]; then
    PRELOAD=true
    if [[ -z "$OUTPUT_STAR" || -z "$STAR_INDEX" ]]; then
        echo "ERROR: output_star_dir and star_index are required when processing FASTQ files"
        exit 1
    fi
fi

echo "Input analysis:"
echo "  FASTQ files detected: $HAS_FASTQ"
echo "  BAM files detected: $HAS_BAM"
echo "  Will preload STAR index: $PRELOAD"

# STAR_INDEX="$HOME/data/gencode/Gencode_human/release_46/STAR_2.7.11b_150"
# STAR_INDEX="$HOME/data/gencode/Gencode_mouse/release_M32/STAR_2.7.10b_150"
# GTF="$HOME/data/gencode/Gencode_mouse/release_M32/gencode.vM32.primary_assembly.annotation.gtf.gz"
# GTF="$HOME/data/gencode/Gencode_human/release_46/gencode.v46.primary_assembly.annotation.gtf.gz"

# Function to create mock featureCounts summary
create_mock_summary() {
    local input_bam="$1"
    local output_summary="$2"
    
    cat > "$output_summary" << EOF
Status	$input_bam
Assigned	0
Unassigned_Unmapped	0
Unassigned_Read_Type	0
Unassigned_Singleton	0
Unassigned_MappingQuality	0
Unassigned_Chimera	0
Unassigned_FragmentLength	0
Unassigned_Duplicate	0
Unassigned_MultiMapping	0
Unassigned_Secondary	0
Unassigned_NonSplit	0
Unassigned_NoFeatures	0
Unassigned_Overlapping_Length	0
Unassigned_Ambiguity	0
EOF
}

# --- Directory Setup ---
FEATURECOUNTS_DIR="$OUTPUT_BAM"
echo "Creating output directories..."
mkdir -p "$FEATURECOUNTS_DIR"

if [[ "$PRELOAD" == true ]]; then
    STAR_DIR="$OUTPUT_STAR"
    mkdir -p "$STAR_DIR"
    
    # --- Pre-load STAR Genome Index ---
    echo "Loading STAR genome into memory..."
    # Remove first?
    STAR \
        --runThreadN 1 \
        --genomeDir "${STAR_INDEX-}" \
        --genomeLoad Remove \
        --outFileNamePrefix _temp/ \
        &>/dev/null || :
    rm -rf _temp 2>/dev/null || :
    STAR \
        --runThreadN 1 \
        --genomeDir "$STAR_INDEX" \
        --genomeLoad LoadAndExit \
        --outFileNamePrefix _temp/ \
        > /dev/null
    rm -rf _temp
fi

# Convert space-separated string back to array
read -ra input_files <<< "$INPUT_FILES"
i=0
for input_file in "${input_files[@]}"; do
    ((++i))
    if [[ $i -le $SKIP ]]; then
        echo "Skipping sample $i: $input_file"
        continue
    fi

    # Determine file type and extract sample name
    if is_alignment_file "$input_file"; then
        # BAM/SAM file processing
        sample_name="$(extract_sample_name "$input_file")"
        echo "Processing BAM/SAM sample: $sample_name"
        
        # For BAM/SAM files, we skip STAR alignment and go directly to featureCounts
        star_bam="$input_file"
        
    else
        # FASTQ file processing
        fqR1="$input_file"
        fqR2="$(derive_r2_from_r1 "$fqR1")"
        sample_name="$(extract_sample_name "$fqR1")"

        # Skip if no matching R2
        if [[ ! -f "$fqR2" ]]; then
            echo "WARNING: Could not find R2 for: $fqR1"
            continue
        fi
        echo "Processing FASTQ sample: $sample_name"

        # Lessons learned during crafting the code below:
        # The big command below basically drops short reads using cutadapt and directly pipes its 
        # output to STAR, without saving any intermediate files.
        # FIFO files are used for communication between processes. To make it work, I had to
        # 1. Explicitly say `exec 3<>"$FIFO_R1"` and `exec 4<>"$FIFO_R2"` after creating the FIFOs. This is 
        #  to ensure that the FIFOs are opened for both reading and writing and to prevent blocking. As a
        #  result, we also need to explicitly close the FIFOs after the cutadapt command (normally the FIFO 
        #  will close when it sees EOF), so that STAR knows that the input is complete. Suprisingly, the way 
        #  to do this is to echo a newline to the FIFOs instead of the moore common `exec 3>&- 4>&-`.
        #  I really don't know why.
        # 2. Run the consumer processes (STAR) first in the background (&) and then the producer (cutadapt).
        #   Otherwise, blocking still occurs, though I don't understand why.
        # 3. Record the PID of the STAR process and wait for it to finish before removing the FIFOs and running
        #   next steps. Otherwise the output of STAR might be incomplete.
        # 4. Somehow cutadapt does not work well with FIFO files. That is, if I set -o and -p to FIFOs, it
        #  will just block. Instead, I let cutadapt output interleaved FASTQ to stdout and used another awk
        #  command to split it back into two FIFOs.
        # 5. FIFO file names must be explicit since STAR reads data based on file name extension. We not only
        #  need the FIFOs to be text (not gzipped), but also need to set --readFilesCommand to cat. Just
        #  setting --readFilesCommand to cat still results in error when the FIFOs are named with .gz.

        # Why don't I just bypass cutadapt and use awk directly to two FIFOs lol? cutadapt 
        # seems to be the source of the many surprises here.

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
            --runThreadN "$((THREADS - 2))" \
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
            --outBAMcompression 1 \
            --outTmpKeep None \
            --quantMode GeneCounts \
            > /dev/null &
        star_pid=$!
        cutadapt \
            -j 1 \
            --interleaved \
            --minimum-length 60:60 \
            "$fqR1" "$fqR2" \
            2> /dev/null \
            | awk -v p1="$FIFO_R1" -v p2="$FIFO_R2" \
               '
                    {
                        # Replace space with underscore only on header lines
                        if ((NR - 1) % 4 == 0) {
                            gsub(" ", "_", $0)
                        }

                        # Write to R1 (first 4 lines) or R2 (next 4 lines) in each 8-line block
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
        star_bam_raw="$STAR_DIR/${sample_name}/Aligned.out.bam"
        star_bam="$STAR_DIR/${sample_name}/${sample_name}.bam"
        mv "$star_bam_raw" "$star_bam"
    fi

    # Common processing for both FASTQ and BAM inputs: featureCounts
    # About featureCounts arguments:
    # -T: Number of threads
    # -p: sequencing data is paired-end
    # -B: Only count reads that are properly paired
    # -C: Do not count read pairs that have two ends mapping to different chromosomes or 
    #   mapping to the same chromosome but on different strands.
    # -s: Strand specificity (0 = unstranded, 1 = stranded, 2 = reversely stranded)
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
    # --countReadPairs: featureCounts would drop R2 alignments, which we DO NOT want.
    # featureCounts input and output order:
    #   1. featureCounts makes sure that input is sorted such that alignments with the same QNAMEs are
    #   adjacent. If the input is already sorted as such, --donotsort can be used to 
    #   skip the sorting step.
    #   2. featureCounts output BAM does not have any orders when there are more than 1 
    #   threads. Alignments with the same QNAMEs can be separated.

    first_bam_line=$(samtools view "$star_bam" | head -n1) || true
    final_bam="$FEATURECOUNTS_DIR/${sample_name}.sorted.bam"
    # featurecounts cannot handle empty input. so check if the bam from star is empty. 
    # If so, just copy it and touch other output files
    if ! echo "$first_bam_line" | grep -q .; then
        echo $first_bam_line >&2
        cp "$star_bam" "$final_bam"
        samtools index "$final_bam"
        echo "Empty bam, featurecounts not run" > "$FEATURECOUNTS_DIR/${sample_name}.log"
        # Mock featurecounts output
        echo "# Program:featureCounts v2.0.3; Command:mock" > "$FEATURECOUNTS_DIR/${sample_name}.txt"
        echo "Geneid\tChr\tStart\tEnd\tStrand\tLength\t$star_bam" >> "$FEATURECOUNTS_DIR/${sample_name}.txt"
        # Mock featurecounts summary
        create_mock_summary "$star_bam" "$FEATURECOUNTS_DIR/${sample_name}.txt.summary"
    else
        FIFO="$FEATURECOUNTS_DIR/$(basename $star_bam).featureCounts.bam"
        if [ -e "$FIFO" ] && [ ! -p "$FIFO" ]; then
            echo "Error: $FIFO exists but is not a FIFO" >&2
            exit 1
        else
            rm -f "$FIFO"
        fi
        mkfifo "$FIFO"
        # If not empty, run featureCounts and sambamba sort
        read SORT_THREADS FC_THREADS < <(python -m folitools.scripts.foli_03_map_utils --total-cores "$((THREADS - 1))")
        
        # Use temporary file for BAM output to allow safe overwriting
        TEMP_FC_TXT="$FEATURECOUNTS_DIR/_${sample_name}.txt"
        TEMP_BAM="$FEATURECOUNTS_DIR/_${sample_name}.sorted.bam"
        FINAL_FC_TXT="$FEATURECOUNTS_DIR/${sample_name}.txt"
        
        # Build featureCounts command with optional overlap flags
        FEATURECOUNTS_CMD="featureCounts \
            -T $((FC_THREADS - 1)) \
            -a $GTF_PATH \
            -o $TEMP_FC_TXT \
            -p \
            -B -C \
            -s $STRAND \
            --donotsort \
            -R BAM \
            --Rpath $FEATURECOUNTS_DIR/"
        
        # Add overlap flags if allow_overlap is true
        if [[ "$ALLOW_OVERLAP" == "True" ]]; then
            FEATURECOUNTS_CMD="$FEATURECOUNTS_CMD -O"
        fi

        # Add multimapping flags if allow_multimapping is true
        if [[ "$ALLOW_MULTIMAPPING" == "True" ]]; then
            FEATURECOUNTS_CMD="$FEATURECOUNTS_CMD -M"
        fi

        # Add fraction flag if either overlap or multimapping is true
        if [[ "$ALLOW_OVERLAP" == "True" || "$ALLOW_MULTIMAPPING" == "True" ]]; then
            FEATURECOUNTS_CMD="$FEATURECOUNTS_CMD --fraction"
        fi
        
        FEATURECOUNTS_CMD="$FEATURECOUNTS_CMD $star_bam"
        
        eval "$FEATURECOUNTS_CMD" \
            2> "$FEATURECOUNTS_DIR/$sample_name.log" \
        & \
        samtools collate \
            -O \
            -l 1 \
            --threads 1 \
            "$FIFO" \
        | \
        python -m folitools.add_tags --cell_tag ${sample_name} \
        | \
        sambamba sort \
            --nthreads "$SORT_THREADS" \
            --memory-limit 16GB \
            --out "$TEMP_BAM" \
            /dev/stdin \
            2> /dev/null \
        & \
        wait

        rm -f "$FIFO"

        # Check if the temporary BAM file was created successfully
        if [[ ! -f "$TEMP_BAM" ]]; then
            echo "ERROR: Failed to create output BAM file: $TEMP_BAM" >&2
            exit 1
        fi
        
        # Move temporary files to final location (allows overwriting)
        mv "$TEMP_BAM" "$final_bam"
        mv "$TEMP_BAM.bai" "$final_bam.bai"
        mv "$TEMP_FC_TXT" "$FINAL_FC_TXT"
        mv "$TEMP_FC_TXT.summary" "$FINAL_FC_TXT.summary"
        
        # Clean up intermediate BAM file for FASTQ inputs (not for BAM/SAM inputs as that's the original)
        if is_fastq_file "$input_file"; then
            rm "$star_bam"
        fi
    fi
    
    # Remove input files to save space
    if [[ "$DELETE" == "True" ]]; then
        if is_alignment_file "$input_file"; then
            # For BAM/SAM inputs, only delete if the input file is different from the final output
            if [[ "$(realpath "$input_file")" != "$(realpath "$final_bam")" ]]; then
                rm "$input_file"
            fi
        else
            rm "$fqR1" "$fqR2"
        fi
    fi

    # - STAR can read from stdin with FIFO, and can output BAM to stdout.
    # - featureCounts cannot read BAM from stdin when -R is specified, but can output BAM with FIFO.
    # - sambamba sort can take input from stdin and output to stdout.
    # - umi_tools group in the next step does not support input from stdin but can 
    #    output BAM to stdout for sorting.
    # - umi_tools group does not require sorted BAM, so I do not sort bam here and 
    #    directly give unsorted BAM to umi_tools group in the next step.
done | tqdm --total ${#input_files[@]} > /dev/null

# Unload the STAR genome index to free up memory (only if we loaded it)
if [[ "$PRELOAD" == true ]]; then
    echo "Unloading STAR genome from memory."
    STAR \
        --runThreadN 1 \
        --genomeDir "$STAR_INDEX" \
        --genomeLoad Remove \
        --outFileNamePrefix _temp/ \
        > /dev/null
    rm -rf _temp
fi

echo "Pipeline finished successfully."