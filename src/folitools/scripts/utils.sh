# File format validation and processing utilities

# Check if a file is a supported FASTQ format
is_fastq_file() {
    local file="$1"
    [[ "$file" == *.fq.gz || "$file" == *.fastq.gz || "$file" == *.fq || "$file" == *.fastq ]]
}

# Check if a file is a supported alignment format  
is_alignment_file() {
    local file="$1"
    [[ "$file" == *.bam || "$file" == *.sam ]]
}

# Extract sample name from any supported file format
extract_sample_name() {
    local file="$1"
    local basename_file sample_name
    
    if [[ "$file" == *.bam ]]; then
        basename_file="$(basename "$file" .bam)"
    elif [[ "$file" == *.sam ]]; then
        basename_file="$(basename "$file" .sam)"
    elif [[ "$file" == *.fq.gz ]]; then
        basename_file="$(basename "$file" .fq.gz)"
    elif [[ "$file" == *.fastq.gz ]]; then
        basename_file="$(basename "$file" .fastq.gz)"
    elif [[ "$file" == *.fq ]]; then
        basename_file="$(basename "$file" .fq)"
    elif [[ "$file" == *.fastq ]]; then
        basename_file="$(basename "$file" .fastq)"
    else
        echo "ERROR: Unsupported file format: $file" >&2
        return 1
    fi
    
    # Extract the sample name (remove R1/R2 indicators if present)
    # First remove common R1/R2 patterns like _R1, _1, .R1, .1, etc.
    sample_name="$basename_file"
    sample_name="${sample_name%_R1}"
    sample_name="${sample_name%_R2}"
    sample_name="${sample_name%_1}"
    sample_name="${sample_name%_2}"
    sample_name="${sample_name%.R1}"
    sample_name="${sample_name%.R2}"
    sample_name="${sample_name%.1}"
    sample_name="${sample_name%.2}"
    echo "$sample_name"
}

# Derive R2 file path from R1 file path (handles multiple patterns)
derive_r2_from_r1() {
    local r1_file="$1"
    local r2_file
    
    # Try different R1/R2 patterns - order matters!
    if [[ "$r1_file" == *_R1.* ]]; then
        r2_file="${r1_file/_R1./_R2.}"
    elif [[ "$r1_file" == *.R1.* ]]; then
        r2_file="${r1_file/.R1./.R2.}"
    elif [[ "$r1_file" == *_1.* ]]; then
        r2_file="${r1_file/_1./_2.}"
    elif [[ "$r1_file" == *.1.* ]]; then
        r2_file="${r1_file/.1./.2.}"
    elif [[ "$r1_file" == *_R1_* ]]; then
        r2_file="${r1_file/_R1_/_R2_}"
    elif [[ "$r1_file" == *_r1.* ]]; then
        r2_file="${r1_file/_r1./_r2.}"
    else
        # If no R1/1 pattern found, try adding _R2 before extension
        local basename="${r1_file%.*}"
        local extension="${r1_file##*.}"
        r2_file="${basename}_R2.${extension}"
    fi
    
    echo "$r2_file"
}

# Validate file formats for a list of files
validate_file_formats() {
    local file_list="$1"
    local allowed_formats="$2"  # "fastq", "alignment", or "both"
    
    read -ra files <<< "$file_list"
    
    for file in "${files[@]}"; do
        case "$allowed_formats" in
            "fastq")
                if ! is_fastq_file "$file"; then
                    echo "ERROR: Unsupported file format: $file. Only FASTQ files (.fq/.fastq/.fq.gz/.fastq.gz) are supported." >&2
                    return 1
                fi
                ;;
            "alignment")
                if ! is_alignment_file "$file"; then
                    echo "ERROR: Unsupported file format: $file. Only alignment files (.bam/.sam) are supported." >&2
                    return 1
                fi
                ;;
            "both")
                if ! is_fastq_file "$file" && ! is_alignment_file "$file"; then
                    echo "ERROR: Unsupported file format: $file. Only FASTQ (.fq/.fastq/.fq.gz/.fastq.gz) and alignment (.bam/.sam) files are supported." >&2
                    return 1
                fi
                ;;
            *)
                echo "ERROR: Invalid allowed_formats parameter: $allowed_formats" >&2
                return 1
                ;;
        esac
    done
    
    return 0
}

run_seqkit_stats() {
    # Usage: run_seqkit_stats /path/to/fastq_dir

    local fastq_dir="$1"
    # convert to absolute path
    fastq_dir=$(realpath "$fastq_dir")
    local stats_file="$fastq_dir.stats"

    if [[ -z "$fastq_dir" ]]; then
        echo "Error: Missing required argument: FASTQ directory path." >&2
        return 1
    fi

    # Find matching FASTQ files
    local fastq_files=()
    mapfile -d '' -t fastq_files < <(
        find "$fastq_dir" -maxdepth 1 -type f \
            \( -name "*.fastq.gz" -o -name "*.fq.gz" -o -name "*.fq" -o -name "*.fastq" \) \
            -print0
    )
    if [[ ${#fastq_files[@]} -eq 0 ]]; then
        echo "Error: No FASTQ files found in '$fastq_dir'." >&2
        return 1
    fi

    # if [[ -f "$stats_file" ]]; then
    #     echo "Warning: Stats file '$stats_file' already exists. Skipping seqkit stats."
    # else
        # sort the files by name to ensure consistent order
        # IFS=$'\n' sorted_files=($(sort <<<"${fastq_files[*]}"))
        # unset IFS
    #     mapfile -d '' -t sorted_files < <(
    #         printf '%s\0' "${fastq_files[@]}" | sort -zV
    #     )
    #     seqkit stats --all --tabular --threads "${THREADS:-1}" "${sorted_files[@]}" > "$stats_file"
    # fi
    mapfile -d '' -t sorted_files < <(
        printf '%s\0' "${fastq_files[@]}" | sort -zV
    )
    seqkit stats --all --tabular --threads "${THREADS:-1}" "${sorted_files[@]}" > "$stats_file"

}


run_fastqc() {
    # Usage: run_fastqc <R1> <R2> <output_dir> <threads>
    # This function effectively just checks if the files are empty and runs fastqc if 
    # they are not, as fastqc will raise error upon empty files.
    local R1="$1"
    local R2="$2"
    local output_dir="$3"
    local threads="$4"

    # Check for required arguments
    if [[ -z "$R1" || -z "$R2" || -z "$output_dir" || -z "$threads" ]]; then
        echo "Error: Missing arguments. Usage: run_fastqc <R1> <R2> <output_dir> <threads>" >&2
        return 1
    fi

    flag=1
    for file in "$R1" "$R2"; do
        if [[ ! -f "$file" ]]; then
            echo "Error: File '$file' does not exist." >&2
            return 1
        fi

        if file "$file" | grep -q 'gzip compressed'; then
            if (( $(gzip -cd -- "$file" | head -c1 | wc -c) == 0 )); then
                flag=0
            fi
        else
            if (( $(head -c1 -- "$file" | wc -c) == 0 )); then
                flag=0
            fi
        fi
    done
    if [[ $flag -eq 0 ]]; then
        echo "Warning: $R1 or $R2 is empty. Skipping FastQC." >&2
        return
    fi

    # Run FastQC
    fastqc -t "$threads" -o "$output_dir" "$R1" "$R2" &> /dev/null
}
