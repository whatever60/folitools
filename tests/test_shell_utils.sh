#!/usr/bin/env bash

# Unit tests for shell utility functions
set -euo pipefail

# Get the script directory and source utils
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
UTILS_DIR="$(dirname "$SCRIPT_DIR")/src/folitools/scripts"
source "$UTILS_DIR/utils.sh"

# Test counters
TESTS_RUN=0
TESTS_PASSED=0
TESTS_FAILED=0

# Test framework functions
assert_equal() {
    local expected="$1"
    local actual="$2"
    local test_name="$3"
    
    ((TESTS_RUN++))
    
    if [[ "$expected" == "$actual" ]]; then
        echo "✓ PASS: $test_name"
        ((TESTS_PASSED++))
    else
        echo "✗ FAIL: $test_name"
        echo "  Expected: '$expected'"
        echo "  Actual:   '$actual'"
        ((TESTS_FAILED++))
    fi
}

assert_true() {
    local test_name="$2"
    
    ((TESTS_RUN++))
    
    # Execute the condition directly instead of using eval
    if $1; then
        echo "✓ PASS: $test_name"
        ((TESTS_PASSED++))
    else
        echo "✗ FAIL: $test_name"
        echo "  Condition '$1' was false"
        ((TESTS_FAILED++))
    fi
}

assert_false() {
    local test_name="$2"
    
    ((TESTS_RUN++))
    
    # Execute the condition directly instead of using eval
    if ! $1; then
        echo "✓ PASS: $test_name"
        ((TESTS_PASSED++))
    else
        echo "✗ FAIL: $test_name"
        echo "  Condition '$1' was true (expected false)"
        ((TESTS_FAILED++))
    fi
}

# Test is_fastq_file function
test_is_fastq_file() {
    echo "Testing is_fastq_file function..."
    
    # Test positive cases
    if is_fastq_file "sample.fq.gz"; then
        echo "✓ PASS: is_fastq_file: .fq.gz"
        ((TESTS_PASSED++))
    else
        echo "✗ FAIL: is_fastq_file: .fq.gz"
        ((TESTS_FAILED++))
    fi
    ((TESTS_RUN++))
    
    if is_fastq_file "sample.fastq.gz"; then
        echo "✓ PASS: is_fastq_file: .fastq.gz"
        ((TESTS_PASSED++))
    else
        echo "✗ FAIL: is_fastq_file: .fastq.gz"
        ((TESTS_FAILED++))
    fi
    ((TESTS_RUN++))
    
    if is_fastq_file "sample.fq"; then
        echo "✓ PASS: is_fastq_file: .fq"
        ((TESTS_PASSED++))
    else
        echo "✗ FAIL: is_fastq_file: .fq"
        ((TESTS_FAILED++))
    fi
    ((TESTS_RUN++))
    
    if is_fastq_file "sample.fastq"; then
        echo "✓ PASS: is_fastq_file: .fastq"
        ((TESTS_PASSED++))
    else
        echo "✗ FAIL: is_fastq_file: .fastq"
        ((TESTS_FAILED++))
    fi
    ((TESTS_RUN++))
    
    # Test negative cases
    if ! is_fastq_file "sample.bam"; then
        echo "✓ PASS: is_fastq_file: .bam (should be false)"
        ((TESTS_PASSED++))
    else
        echo "✗ FAIL: is_fastq_file: .bam (should be false)"
        ((TESTS_FAILED++))
    fi
    ((TESTS_RUN++))
    
    if ! is_fastq_file "sample.sam"; then
        echo "✓ PASS: is_fastq_file: .sam (should be false)"
        ((TESTS_PASSED++))
    else
        echo "✗ FAIL: is_fastq_file: .sam (should be false)"
        ((TESTS_FAILED++))
    fi
    ((TESTS_RUN++))
    
    if ! is_fastq_file "sample.txt"; then
        echo "✓ PASS: is_fastq_file: .txt (should be false)"
        ((TESTS_PASSED++))
    else
        echo "✗ FAIL: is_fastq_file: .txt (should be false)"
        ((TESTS_FAILED++))
    fi
    ((TESTS_RUN++))
}

# Test is_alignment_file function
test_is_alignment_file() {
    echo "Testing is_alignment_file function..."
    
    assert_true 'is_alignment_file "sample.bam"' "is_alignment_file: .bam"
    assert_true 'is_alignment_file "sample.sam"' "is_alignment_file: .sam"
    
    assert_false 'is_alignment_file "sample.fq.gz"' "is_alignment_file: .fq.gz (should be false)"
    assert_false 'is_alignment_file "sample.fastq"' "is_alignment_file: .fastq (should be false)"
    assert_false 'is_alignment_file "sample.txt"' "is_alignment_file: .txt (should be false)"
}

# Test extract_sample_name function
test_extract_sample_name() {
    echo "Testing extract_sample_name function..."
    
    # FASTQ files
    assert_equal "sample1" "$(extract_sample_name "sample1_1.fq.gz")" "extract_sample_name: sample1_1.fq.gz"
    assert_equal "sample1" "$(extract_sample_name "sample1_R1.fastq.gz")" "extract_sample_name: sample1_R1.fastq.gz"
    assert_equal "sample1" "$(extract_sample_name "sample1.fq")" "extract_sample_name: sample1.fq"
    assert_equal "sample1" "$(extract_sample_name "sample1.fastq")" "extract_sample_name: sample1.fastq"
    
    # Alignment files
    assert_equal "sample1" "$(extract_sample_name "sample1.sorted.bam")" "extract_sample_name: sample1.sorted.bam"
    assert_equal "sample1" "$(extract_sample_name "sample1.bam")" "extract_sample_name: sample1.bam"
    assert_equal "sample1" "$(extract_sample_name "sample1.sam")" "extract_sample_name: sample1.sam"
    
    # Complex sample names
    assert_equal "CTRL-01" "$(extract_sample_name "CTRL-01_1.fq.gz")" "extract_sample_name: CTRL-01_1.fq.gz"
    assert_equal "sample" "$(extract_sample_name "sample.with.dots_1.fq.gz")" "extract_sample_name: sample.with.dots_1.fq.gz"
}

# Test derive_r2_from_r1 function
test_derive_r2_from_r1() {
    echo "Testing derive_r2_from_r1 function..."
    
    # Standard _1/_2 pattern
    assert_equal "sample_2.fq.gz" "$(derive_r2_from_r1 "sample_1.fq.gz")" "derive_r2_from_r1: _1/_2 pattern"
    
    # _R1_/_R2_ pattern (fastp style)
    assert_equal "sample_R2_.fastq.gz" "$(derive_r2_from_r1 "sample_R1_.fastq.gz")" "derive_r2_from_r1: _R1_/_R2_ pattern"
    
    # _r1./_r2. pattern
    assert_equal "sample_r2.fq" "$(derive_r2_from_r1 "sample_r1.fq")" "derive_r2_from_r1: _r1./_r2. pattern"
    
    # _R1./_R2. pattern
    assert_equal "sample_R2.fastq" "$(derive_r2_from_r1 "sample_R1.fastq")" "derive_r2_from_r1: _R1./_R2. pattern"
    
    # Fallback to _1/_2 when pattern not recognized
    assert_equal "sample_other_2.fq.gz" "$(derive_r2_from_r1 "sample_other_1.fq.gz")" "derive_r2_from_r1: fallback pattern"
}

# Test validate_file_formats function
test_validate_file_formats() {
    echo "Testing validate_file_formats function..."
    
    # Test fastq validation
    assert_true 'validate_file_formats "sample1.fq.gz sample2.fastq" "fastq"' "validate_file_formats: valid fastq files"
    
    # Test alignment validation  
    assert_true 'validate_file_formats "sample1.bam sample2.sam" "alignment"' "validate_file_formats: valid alignment files"
    
    # Test both validation
    assert_true 'validate_file_formats "sample1.fq.gz sample2.bam" "both"' "validate_file_formats: mixed valid files"
    
    # Test validation failures (these should return false, but we can't easily test stderr)
    # Note: These tests will produce error messages, which is expected
    echo "  Note: The following tests will show expected error messages..."
    
    # Redirect stderr to suppress expected error messages in our test output
    assert_false 'validate_file_formats "sample1.txt" "fastq" 2>/dev/null' "validate_file_formats: invalid fastq file"
    assert_false 'validate_file_formats "sample1.fq.gz" "alignment" 2>/dev/null' "validate_file_formats: fastq file with alignment validation"
}

# Test integration with actual file processing scenarios
test_integration_scenarios() {
    echo "Testing integration scenarios..."
    
    # Create temporary test files
    local temp_dir="/tmp/shell_utils_test_$$"
    mkdir -p "$temp_dir"
    
    # Create test files with different naming patterns
    touch "$temp_dir/sample1_1.fq.gz"
    touch "$temp_dir/sample1_2.fq.gz"
    touch "$temp_dir/sample2_R1_.fastq.gz"
    touch "$temp_dir/sample2_R2_.fastq.gz"
    touch "$temp_dir/sample3.bam"
    touch "$temp_dir/sample4.sam"
    
    # Test sample name extraction from real files
    assert_equal "sample1" "$(extract_sample_name "$temp_dir/sample1_1.fq.gz")" "integration: sample1 extraction"
    assert_equal "sample2" "$(extract_sample_name "$temp_dir/sample2_R1_.fastq.gz")" "integration: sample2 extraction"
    assert_equal "sample3" "$(extract_sample_name "$temp_dir/sample3.bam")" "integration: sample3 extraction"
    assert_equal "sample4" "$(extract_sample_name "$temp_dir/sample4.sam")" "integration: sample4 extraction"
    
    # Test R2 derivation with file existence check
    local r2_file="$(derive_r2_from_r1 "$temp_dir/sample1_1.fq.gz")"
    assert_equal "$temp_dir/sample1_2.fq.gz" "$r2_file" "integration: R2 derivation"
    assert_true '[[ -f "$r2_file" ]]' "integration: derived R2 file exists"
    
    # Clean up
    rm -rf "$temp_dir"
}

# Run all tests
main() {
    echo "Running shell utility tests..."
    echo "============================="
    
    test_is_fastq_file
    echo
    test_is_alignment_file  
    echo
    test_extract_sample_name
    echo
    test_derive_r2_from_r1
    echo
    test_validate_file_formats
    echo
    test_integration_scenarios
    
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
}

# Run tests if script is executed directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi
