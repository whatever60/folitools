#!/usr/bin/env bash

# Simple unit tests for shell utility functions
set -euo pipefail

# Get the script directory and source utils
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
UTILS_DIR="$(dirname "$SCRIPT_DIR")/src/folitools/scripts"
source "$UTILS_DIR/utils.sh"

# Test counters
TESTS_RUN=0
TESTS_PASSED=0
TESTS_FAILED=0

# Simple test function
run_test() {
    local test_name="$1"
    local command="$2"
    
    ((TESTS_RUN++))
    
    if eval "$command"; then
        echo "✓ PASS: $test_name"
        ((TESTS_PASSED++))
    else
        echo "✗ FAIL: $test_name"
        ((TESTS_FAILED++))
    fi
}

echo "Running shell utility tests..."
echo "============================="

# Test is_fastq_file function
echo "Testing is_fastq_file function..."
run_test "is_fastq_file: .fq.gz" 'is_fastq_file "sample.fq.gz"'
run_test "is_fastq_file: .fastq.gz" 'is_fastq_file "sample.fastq.gz"'
run_test "is_fastq_file: .fq" 'is_fastq_file "sample.fq"'
run_test "is_fastq_file: .fastq" 'is_fastq_file "sample.fastq"'
run_test "is_fastq_file: .bam (should fail)" '! is_fastq_file "sample.bam"'
run_test "is_fastq_file: .sam (should fail)" '! is_fastq_file "sample.sam"'

echo
# Test is_alignment_file function  
echo "Testing is_alignment_file function..."
run_test "is_alignment_file: .bam" 'is_alignment_file "sample.bam"'
run_test "is_alignment_file: .sam" 'is_alignment_file "sample.sam"'
run_test "is_alignment_file: .fq.gz (should fail)" '! is_alignment_file "sample.fq.gz"'

echo
# Test extract_sample_name function
echo "Testing extract_sample_name function..."
run_test "extract_sample_name: sample1_1.fq.gz" '[[ "$(extract_sample_name "sample1_1.fq.gz")" == "sample1" ]]'
run_test "extract_sample_name: sample1.bam" '[[ "$(extract_sample_name "sample1.bam")" == "sample1" ]]'
run_test "extract_sample_name: sample1.sam" '[[ "$(extract_sample_name "sample1.sam")" == "sample1" ]]'
run_test "extract_sample_name: CTRL-01_1.fq.gz" '[[ "$(extract_sample_name "CTRL-01_1.fq.gz")" == "CTRL-01" ]]'

echo  
# Test derive_r2_from_r1 function
echo "Testing derive_r2_from_r1 function..."
run_test "derive_r2_from_r1: _1/_2 pattern" '[[ "$(derive_r2_from_r1 "sample_1.fq.gz")" == "sample_2.fq.gz" ]]'
run_test "derive_r2_from_r1: _R1_/_R2_ pattern" '[[ "$(derive_r2_from_r1 "sample_R1_.fastq.gz")" == "sample_R2_.fastq.gz" ]]'

echo
# Test validate_file_formats function
echo "Testing validate_file_formats function..."
run_test "validate_file_formats: valid fastq files" 'validate_file_formats "sample1.fq.gz sample2.fastq" "fastq" 2>/dev/null'
run_test "validate_file_formats: valid alignment files" 'validate_file_formats "sample1.bam sample2.sam" "alignment" 2>/dev/null'
run_test "validate_file_formats: mixed valid files" 'validate_file_formats "sample1.fq.gz sample2.bam" "both" 2>/dev/null'
run_test "validate_file_formats: invalid file (should fail)" '! validate_file_formats "sample1.txt" "fastq" 2>/dev/null'

echo
echo "============================="
echo "Test Results:"
echo "  Total: $TESTS_RUN"
echo "  Passed: $TESTS_PASSED"
echo "  Failed: $TESTS_FAILED"

if [[ $TESTS_FAILED -eq 0 ]]; then
    echo "✓ All tests passed!"
    exit 0
else
    echo "✗ Some tests failed!"
    exit 1
fi
