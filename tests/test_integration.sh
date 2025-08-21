#!/usr/bin/env bash

# Integration tests for the improved shell scripts
set -euo pipefail

# Test directory
TEST_DIR="/tmp/folitools_integration_test_$$"
mkdir -p "$TEST_DIR"

# Cleanup function
cleanup() {
    rm -rf "$TEST_DIR"
}
trap cleanup EXIT

# Test counters
TESTS_RUN=0
TESTS_PASSED=0
TESTS_FAILED=0

# Test function
test_result() {
    local test_name="$1"
    local success="$2"
    
    ((TESTS_RUN++))
    
    if [[ "$success" == "true" ]]; then
        echo "✓ PASS: $test_name"
        ((TESTS_PASSED++))
    else
        echo "✗ FAIL: $test_name"
        ((TESTS_FAILED++))
    fi
}

echo "Running integration tests for improved shell scripts..."
echo "====================================================="

# Create test files with different formats and naming patterns
echo "Setting up test files..."

# Different FASTQ formats
echo "sample data" > "$TEST_DIR/sample1_1.fq.gz"
echo "sample data" > "$TEST_DIR/sample1_2.fq.gz"
echo "sample data" > "$TEST_DIR/sample2_R1_.fastq.gz"
echo "sample data" > "$TEST_DIR/sample2_R2_.fastq.gz"
echo "sample data" > "$TEST_DIR/sample3.fq"
echo "sample data" > "$TEST_DIR/sample3_2.fq"
echo "sample data" > "$TEST_DIR/sample4.fastq"
echo "sample data" > "$TEST_DIR/sample4_2.fastq"

# Alignment files
echo "sample bam data" > "$TEST_DIR/sample5.bam"
echo "sample sam data" > "$TEST_DIR/sample6.sam"
echo "sample bam data" > "$TEST_DIR/sample7.sorted.bam"

# Create a dummy GTF file
echo "dummy gtf" > "$TEST_DIR/dummy.gtf"

echo "Created test files:"
ls -la "$TEST_DIR"
echo

# Test each script's sample name extraction by checking their help/validation
echo "Testing script validation and help messages..."

# Test foli_01_fastp.sh - should only accept FASTQ files
echo "Testing foli_01_fastp.sh validation..."
if bash /home/ubuntu/dev/folitools/src/folitools/scripts/foli_01_fastp.sh "$TEST_DIR/sample1_1.fq.gz" "$TEST_DIR/out1" 1 0 false 2>&1 | grep -q "ERROR.*Unsupported file format"; then
    test_result "foli_01_fastp.sh: rejects non-FASTQ files" "false"
else
    # This script should accept the FASTQ file (but might fail later for other reasons)
    test_result "foli_01_fastp.sh: accepts FASTQ files" "true"
fi

# Try with a BAM file - should be rejected
if bash /home/ubuntu/dev/folitools/src/folitools/scripts/foli_01_fastp.sh "$TEST_DIR/sample5.bam" "$TEST_DIR/out1" 1 0 false 2>&1 | grep -q "ERROR.*Unsupported file format"; then
    test_result "foli_01_fastp.sh: correctly rejects BAM files" "true"
else
    test_result "foli_01_fastp.sh: correctly rejects BAM files" "false"
fi

# Test foli_02_cutadapt.sh - should only accept FASTQ files
echo
echo "Testing foli_02_cutadapt.sh validation..."
if bash /home/ubuntu/dev/folitools/src/folitools/scripts/foli_02_cutadapt.sh "$TEST_DIR/sample5.bam" "$TEST_DIR/out2" "$TEST_DIR/dummy.gtf" "$TEST_DIR/dummy.gtf" 1 0 false 2>&1 | grep -q "ERROR.*Unsupported file format"; then
    test_result "foli_02_cutadapt.sh: correctly rejects BAM files" "true"
else
    test_result "foli_02_cutadapt.sh: correctly rejects BAM files" "false"
fi

# Test foli_04_count.sh - should only accept alignment files
echo
echo "Testing foli_04_count.sh validation..."
if bash /home/ubuntu/dev/folitools/src/folitools/scripts/foli_04_count.sh "$TEST_DIR/sample1_1.fq.gz" "$TEST_DIR/out4" 1 0 2>&1 | grep -q "ERROR.*Unsupported file format"; then
    test_result "foli_04_count.sh: correctly rejects FASTQ files" "true"
else
    test_result "foli_04_count.sh: correctly rejects FASTQ files" "false"
fi

# Test foli_03_map.sh - should accept both FASTQ and alignment files
echo
echo "Testing foli_03_map.sh validation..."

# Test with mixed files - should accept them but require STAR params for FASTQ
if bash /home/ubuntu/dev/folitools/src/folitools/scripts/foli_03_map.sh "$TEST_DIR/sample1_1.fq.gz $TEST_DIR/sample5.bam" "$TEST_DIR/out3_bam" "" "" "$TEST_DIR/dummy.gtf" 1 0 false 0 2>&1 | grep -q "ERROR.*output_star_dir and star_index are required"; then
    test_result "foli_03_map.sh: correctly requires STAR params for FASTQ input" "true"
else
    test_result "foli_03_map.sh: correctly requires STAR params for FASTQ input" "false"
fi

# Test with only BAM files - should not require STAR params
if bash /home/ubuntu/dev/folitools/src/folitools/scripts/foli_03_map.sh "$TEST_DIR/sample5.bam" "$TEST_DIR/out3_bam" "" "" "$TEST_DIR/dummy.gtf" 1 0 false 0 2>&1 | grep -q "Input analysis"; then
    test_result "foli_03_map.sh: accepts BAM-only input without STAR params" "true"
else
    test_result "foli_03_map.sh: accepts BAM-only input without STAR params" "false"
fi

# Test unsupported file format
if bash /home/ubuntu/dev/folitools/src/folitools/scripts/foli_03_map.sh "$TEST_DIR/dummy.gtf" "$TEST_DIR/out3_bam" "" "" "$TEST_DIR/dummy.gtf" 1 0 false 0 2>&1 | grep -q "ERROR.*Unsupported file format"; then
    test_result "foli_03_map.sh: correctly rejects unsupported file formats" "true"
else
    test_result "foli_03_map.sh: correctly rejects unsupported file formats" "false"
fi

echo
echo "Testing utility functions through scripts..."

# Test that different file extensions are handled properly by sourcing utils directly
cd /home/ubuntu/dev/folitools
source src/folitools/scripts/utils.sh

# Test sample name extraction
sample_name_fq=$(extract_sample_name "$TEST_DIR/sample1_1.fq.gz" 2>/dev/null || echo "FAILED")
if [[ "$sample_name_fq" == "sample1" ]]; then
    test_result "extract_sample_name: sample1_1.fq.gz -> sample1" "true"
else
    test_result "extract_sample_name: sample1_1.fq.gz -> sample1" "false"
fi

sample_name_bam=$(extract_sample_name "$TEST_DIR/sample7.sorted.bam" 2>/dev/null || echo "FAILED")
if [[ "$sample_name_bam" == "sample7" ]]; then
    test_result "extract_sample_name: sample7.sorted.bam -> sample7" "true"
else
    test_result "extract_sample_name: sample7.sorted.bam -> sample7" "false"
fi

# Test R2 derivation
r2_file=$(derive_r2_from_r1 "$TEST_DIR/sample1_1.fq.gz" 2>/dev/null || echo "FAILED")
if [[ "$r2_file" == "$TEST_DIR/sample1_2.fq.gz" ]]; then
    test_result "derive_r2_from_r1: _1/_2 pattern works" "true"
else
    test_result "derive_r2_from_r1: _1/_2 pattern works" "false"
fi

r2_file_r=$(derive_r2_from_r1 "$TEST_DIR/sample2_R1_.fastq.gz" 2>/dev/null || echo "FAILED")
if [[ "$r2_file_r" == "$TEST_DIR/sample2_R2_.fastq.gz" ]]; then
    test_result "derive_r2_from_r1: _R1_/_R2_ pattern works" "true"
else
    test_result "derive_r2_from_r1: _R1_/_R2_ pattern works" "false"
fi

# Test file type detection
if is_fastq_file "$TEST_DIR/sample1_1.fq.gz" 2>/dev/null; then
    test_result "is_fastq_file: correctly identifies .fq.gz" "true"
else
    test_result "is_fastq_file: correctly identifies .fq.gz" "false"
fi

if is_alignment_file "$TEST_DIR/sample5.bam" 2>/dev/null; then
    test_result "is_alignment_file: correctly identifies .bam" "true"
else
    test_result "is_alignment_file: correctly identifies .bam" "false"
fi

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
