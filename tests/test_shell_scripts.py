"""
Unified tests for shell scripts using pytest.

This module provides comprehensive testing for all shell scripts in the folitools package,
integrating them into the main pytest framework for unified test running and reporting.
"""

import subprocess
import tempfile
from pathlib import Path
from typing import List, Tuple

import pytest


# Test fixtures
@pytest.fixture
def scripts_dir():
    """Get the path to the shell scripts directory."""
    return Path(__file__).parent.parent / "src" / "folitools" / "scripts"


@pytest.fixture
def temp_test_dir():
    """Create a temporary directory for test files."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        yield Path(tmp_dir)


@pytest.fixture
def sample_test_files(temp_test_dir):
    """Create sample test files for validation."""
    files = {
        "fastq": [
            "sample_001_R1.fq.gz",
            "sample_002_R1.fastq.gz",
            "patient_A.1.fq",
            "ctrl.R1.fastq",
        ],
        "alignment": ["sample_001.bam", "sample_002.sam", "result.bam"],
    }

    # Create empty files for testing
    created_files = {"fastq": [], "alignment": []}
    for file_type, file_list in files.items():
        for filename in file_list:
            file_path = temp_test_dir / filename
            file_path.touch()
            created_files[file_type].append(str(file_path))

    return created_files


class TestUtilityFunctions:
    """Test the utility functions in utils.sh."""

    def run_bash_function(
        self, scripts_dir: Path, function_name: str, args: List[str]
    ) -> Tuple[int, str, str]:
        """Helper to run a bash function and return exit code, stdout, stderr."""
        cmd = [
            "bash",
            "-c",
            f"source {scripts_dir}/utils.sh && {function_name} {' '.join(args)}",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        return result.returncode, result.stdout.strip(), result.stderr.strip()

    def test_is_fastq_file(self, scripts_dir):
        """Test FASTQ file detection."""
        test_cases = [
            ("test.fq", 0),
            ("test.fastq", 0),
            ("test.fq.gz", 0),
            ("test.fastq.gz", 0),
            ("test.bam", 1),
            ("test.sam", 1),
            ("test.txt", 1),
        ]

        for filename, expected_exit_code in test_cases:
            exit_code, stdout, stderr = self.run_bash_function(
                scripts_dir, "is_fastq_file", [filename]
            )
            assert exit_code == expected_exit_code, (
                f"Failed for {filename}: expected {expected_exit_code}, got {exit_code}"
            )

    def test_is_alignment_file(self, scripts_dir):
        """Test alignment file detection."""
        test_cases = [
            ("test.bam", 0),
            ("test.sam", 0),
            ("test.fq", 1),
            ("test.fastq", 1),
            ("test.fq.gz", 1),
            ("test.txt", 1),
        ]

        for filename, expected_exit_code in test_cases:
            exit_code, stdout, stderr = self.run_bash_function(
                scripts_dir, "is_alignment_file", [filename]
            )
            assert exit_code == expected_exit_code, (
                f"Failed for {filename}: expected {expected_exit_code}, got {exit_code}"
            )

    @pytest.mark.parametrize(
        "input_file,expected_sample",
        [
            ("ABC001_R1.fq.gz", "ABC001"),
            ("sample_123_R2.fastq", "sample_123"),
            ("patient-45.1.fq.gz", "patient-45"),
            ("ctrl.R1.fastq.gz", "ctrl"),
            ("treatment_group_001.fq", "treatment_group_001"),
            ("file_name_with_underscores_R1.fastq.gz", "file_name_with_underscores"),
            ("simple.bam", "simple"),
            ("complex_sample_name.sam", "complex_sample_name"),
        ],
    )
    def test_extract_sample_name(self, scripts_dir, input_file, expected_sample):
        """Test sample name extraction with various file patterns."""
        exit_code, stdout, stderr = self.run_bash_function(
            scripts_dir, "extract_sample_name", [input_file]
        )
        assert exit_code == 0, f"Function failed: {stderr}"
        assert stdout == expected_sample, (
            f"Expected '{expected_sample}', got '{stdout}'"
        )

    @pytest.mark.parametrize(
        "r1_file,expected_r2",
        [
            ("sample_001_R1.fq.gz", "sample_001_R2.fq.gz"),
            ("sample_002_1.fastq.gz", "sample_002_2.fastq.gz"),
            ("sample_003.R1.fq", "sample_003.R2.fq"),
            ("sample_004.1.fastq", "sample_004.2.fastq"),
            (
                "file_name_with_underscores_R1.fastq.gz",
                "file_name_with_underscores_R2.fastq.gz",
            ),
        ],
    )
    def test_derive_r2_from_r1(self, scripts_dir, r1_file, expected_r2):
        """Test R2 file derivation from R1 files."""
        exit_code, stdout, stderr = self.run_bash_function(
            scripts_dir, "derive_r2_from_r1", [r1_file]
        )
        assert exit_code == 0, f"Function failed: {stderr}"
        assert stdout == expected_r2, f"Expected '{expected_r2}', got '{stdout}'"

    def test_validate_file_formats(self, scripts_dir):
        """Test file format validation function."""
        # Test valid FASTQ files
        exit_code, stdout, stderr = self.run_bash_function(
            scripts_dir,
            "validate_file_formats",
            ['"test1.fq.gz test2.fastq test3.fq"', '"fastq"'],
        )
        assert exit_code == 0, f"FASTQ validation failed: {stderr}"

        # Test valid alignment files
        exit_code, stdout, stderr = self.run_bash_function(
            scripts_dir,
            "validate_file_formats",
            ['"test1.bam test2.sam"', '"alignment"'],
        )
        assert exit_code == 0, f"Alignment validation failed: {stderr}"

        # Test mixed files (should pass with "both")
        exit_code, stdout, stderr = self.run_bash_function(
            scripts_dir,
            "validate_file_formats",
            ['"test1.fq.gz test2.bam test3.fastq"', '"both"'],
        )
        assert exit_code == 0, f"Mixed validation failed: {stderr}"

        # Test invalid format (should fail)
        exit_code, stdout, stderr = self.run_bash_function(
            scripts_dir, "validate_file_formats", ['"test1.txt test2.fq.gz"', '"fastq"']
        )
        assert exit_code == 1, "Invalid format validation should have failed"


class TestShellScriptSyntax:
    """Test that all shell scripts have valid syntax."""

    @pytest.mark.parametrize(
        "script_name",
        [
            "utils.sh",
            "foli_01_fastp.sh",
            "foli_02_cutadapt.sh",
            "foli_03_map.sh",
            "foli_04_count.sh",
        ],
    )
    def test_script_syntax(self, scripts_dir, script_name):
        """Test that each script has valid bash syntax."""
        script_path = scripts_dir / script_name
        assert script_path.exists(), f"Script {script_name} not found"

        result = subprocess.run(
            ["bash", "-n", str(script_path)], capture_output=True, text=True
        )
        assert result.returncode == 0, f"Syntax error in {script_name}: {result.stderr}"


class TestShellScriptHelp:
    """Test help functionality in shell scripts."""

    @pytest.mark.parametrize(
        "script_name",
        [
            "foli_01_fastp.sh",
            "foli_02_cutadapt.sh",
            "foli_03_map.sh",
            "foli_04_count.sh",
        ],
    )
    def test_help_function(self, scripts_dir, script_name):
        """Test that help functions work correctly."""
        script_path = scripts_dir / script_name

        # Test -h flag
        result = subprocess.run(
            ["bash", str(script_path), "-h"],
            capture_output=True,
            text=True,
            cwd=scripts_dir,
        )
        assert result.returncode == 0, (
            f"Help function failed for {script_name}: {result.stderr}"
        )
        assert "Usage:" in result.stdout, f"Help output missing for {script_name}"

        # Test --help flag
        result = subprocess.run(
            ["bash", str(script_path), "--help"],
            capture_output=True,
            text=True,
            cwd=scripts_dir,
        )
        assert result.returncode == 0, (
            f"Help function failed for {script_name}: {result.stderr}"
        )
        assert "Usage:" in result.stdout, f"Help output missing for {script_name}"


class TestShellScriptIntegration:
    """Integration tests for shell scripts with actual file operations."""

    def test_fastp_file_validation(self, scripts_dir, sample_test_files, temp_test_dir):
        """Test that foli_01_fastp.sh properly validates input files."""
        # Test with valid FASTQ files
        fastq_files = " ".join(sample_test_files["fastq"])
        output_dir = temp_test_dir / "fastp_output"

        result = subprocess.run(
            [
                "bash",
                str(scripts_dir / "foli_01_fastp.sh"),
                fastq_files,
                str(output_dir),
            ],
            capture_output=True,
            text=True,
            cwd=scripts_dir,
        )
        # Should fail because fastp isn't installed, but file validation should pass
        # Look for specific error messages to ensure it got past validation
        assert "ERROR: Unsupported file format" not in result.stderr

        # Test with invalid files (should fail validation)
        invalid_files = " ".join(sample_test_files["alignment"])
        result = subprocess.run(
            [
                "bash",
                str(scripts_dir / "foli_01_fastp.sh"),
                invalid_files,
                str(output_dir),
            ],
            capture_output=True,
            text=True,
            cwd=scripts_dir,
        )
        assert result.returncode != 0
        assert "ERROR: Unsupported file format" in result.stderr

    def test_count_file_validation(self, scripts_dir, sample_test_files, temp_test_dir):
        """Test that foli_04_count.sh properly validates input files."""
        # Test with valid alignment files
        alignment_files = " ".join(sample_test_files["alignment"])
        output_dir = temp_test_dir / "count_output"

        result = subprocess.run(
            [
                "bash",
                str(scripts_dir / "foli_04_count.sh"),
                alignment_files,
                str(output_dir),
            ],
            capture_output=True,
            text=True,
            cwd=scripts_dir,
        )
        # Should fail because umi_tools isn't installed, but file validation should pass
        assert "ERROR: Unsupported file format" not in result.stderr

        # Test with invalid files (should fail validation)
        invalid_files = " ".join(sample_test_files["fastq"])
        result = subprocess.run(
            [
                "bash",
                str(scripts_dir / "foli_04_count.sh"),
                invalid_files,
                str(output_dir),
            ],
            capture_output=True,
            text=True,
            cwd=scripts_dir,
        )
        assert result.returncode != 0
        assert "ERROR: Unsupported file format" in result.stderr


class TestShellScriptRobustness:
    """Test edge cases and robustness of shell scripts."""

    def test_sample_name_extraction_edge_cases(self, scripts_dir):
        """Test sample name extraction with edge cases."""
        edge_cases = [
            # Complex sample names
            (
                "very_long_sample_name_with_many_underscores_R1.fastq.gz",
                "very_long_sample_name_with_many_underscores",
            ),
            ("SAMPLE-WITH-DASHES_R2.fq", "SAMPLE-WITH-DASHES"),
            ("123numeric_start_R1.fq.gz", "123numeric_start"),
            ("CamelCase_R1.fastq", "CamelCase"),
            # Files without R1/R2 indicators
            ("just_sample_name.fq.gz", "just_sample_name"),
            ("another.sample.bam", "another.sample"),
            # Multiple dots and underscores
            ("sample.v1.2_final_R1.fq.gz", "sample.v1.2_final"),
        ]

        for input_file, expected_sample in edge_cases:
            exit_code, stdout, stderr = self.run_bash_function(
                scripts_dir, "extract_sample_name", [input_file]
            )
            assert exit_code == 0, f"Function failed for {input_file}: {stderr}"
            assert stdout == expected_sample, (
                f"Expected '{expected_sample}', got '{stdout}' for {input_file}"
            )

    def run_bash_function(
        self, scripts_dir: Path, function_name: str, args: List[str]
    ) -> Tuple[int, str, str]:
        """Helper to run a bash function and return exit code, stdout, stderr."""
        cmd = [
            "bash",
            "-c",
            f"source {scripts_dir}/utils.sh && {function_name} {' '.join(args)}",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        return result.returncode, result.stdout.strip(), result.stderr.strip()


@pytest.mark.slow
class TestShellScriptPerformance:
    """Performance and stress tests for shell scripts (marked as slow)."""

    def test_large_file_list_handling(self, scripts_dir, temp_test_dir):
        """Test handling of large numbers of files."""
        # Create many test files
        file_list = []
        for i in range(100):
            filename = f"sample_{i:03d}_R1.fq.gz"
            file_path = temp_test_dir / filename
            file_path.touch()
            file_list.append(str(file_path))

        # Test file format validation with large list
        file_string = " ".join(file_list)
        exit_code, stdout, stderr = self.run_bash_function(
            scripts_dir, "validate_file_formats", [f'"{file_string}"', '"fastq"']
        )
        assert exit_code == 0, f"Large file list validation failed: {stderr}"

    def run_bash_function(
        self, scripts_dir: Path, function_name: str, args: List[str]
    ) -> Tuple[int, str, str]:
        """Helper to run a bash function and return exit code, stdout, stderr."""
        cmd = [
            "bash",
            "-c",
            f"source {scripts_dir}/utils.sh && {function_name} {' '.join(args)}",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        return result.returncode, result.stdout.strip(), result.stderr.strip()


if __name__ == "__main__":
    # Allow running this file directly for quick testing
    pytest.main([__file__, "-v"])
