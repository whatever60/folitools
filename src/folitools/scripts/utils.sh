run_seqkit_stats() {
    # Usage: run_seqkit_stats /path/to/fastq_dir

    local fastq_dir="$1"
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

    if [[ -f "$stats_file" ]]; then
        echo "Warning: Stats file '$stats_file' already exists. Skipping seqkit stats."
    else
        # sort the files by name to ensure consistent order
        IFS=$'\n' sorted_files=($(sort <<<"${fastq_files[*]}"))
        unset IFS
        seqkit stats --all --tabular --threads "${THREADS:-1}" "${sorted_files[@]}" > "$stats_file"
    fi
}


run_fastqc() {
    # Usage: run_fastqc <trimmed_R1> <trimmed_R2> <output_dir> <threads>
    # This function effectively just checks if the files are empty and runs fastqc if 
    # they are not, as fastqc will raise error upon empty files.
    local trimmed_R1="$1"
    local trimmed_R2="$2"
    local output_dir="$3"
    local threads="$4"

    # Check for required arguments
    if [[ -z "$trimmed_R1" || -z "$trimmed_R2" || -z "$output_dir" || -z "$threads" ]]; then
        echo "Error: Missing arguments. Usage: run_fastqc <trimmed_R1> <trimmed_R2> <output_dir> <threads>" >&2
        return 1
    fi

    flag=1
    for file in "$trimmed_R1" "$trimmed_R2"; do
        if [[ ! -f "$file" ]]; then
            echo "Error: File '$file' does not exist." >&2
            return 1
        fi

        if file "$file" | grep -q 'gzip compressed'; then
            if zcat "$file" | head -c1 | grep -q .; then
                flag=0
            fi
        else
            if head -c1 "$file" | grep -q .; then
                flag=0
            fi
        fi
    done

    if [[ $flag -eq 1 ]]; then
        echo "Warning: File '$file' is empty. Skipping FastQC." >&2
        return
    fi

    # Run FastQC
    fastqc -t "$threads" -o "$output_dir" "$trimmed_R1" "$trimmed_R2" &> /dev/null
}
