#!/usr/bin/env python3
"""
Simple test for recover CLI command with IDT file.
"""

import pytest
import pandas as pd
import tempfile
import shutil
from pathlib import Path


class TestRecoverCLI:
    """Test recover CLI command with IDT file."""

    @pytest.fixture
    def idt_test_file(self):
        """Use the small IDT test file."""
        idt_file = Path(__file__).parent / "test_idt_order.xlsx"
        if not idt_file.exists():
            pytest.skip(f"IDT test file not found: {idt_file}")
        return idt_file

    @pytest.fixture
    def temp_output_dir(self):
        """Create temporary output directory."""
        temp_dir = Path(tempfile.mkdtemp(prefix="test_recover_cli_"))
        yield temp_dir
        # Cleanup
        if temp_dir.exists():
            shutil.rmtree(temp_dir)

    def test_recover_cli_help(self):
        """Test that recover CLI shows help without errors."""
        import subprocess
        import sys

        # Test the CLI help
        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "folitools.primer_selection.cli",
                "recover",
                "--help",
            ],
            capture_output=True,
            text=True,
            cwd="/home/ubuntu/dev/folitools",
        )

        # Help should exit with code 0
        assert result.returncode == 0, (
            f"Help should exit with code 0, got {result.returncode}"
        )

        help_text = result.stdout
        assert "recover" in help_text.lower(), "Help should mention recover command"
        assert "order-excel" in help_text.lower(), (
            "Help should mention order-excel parameter"
        )

    def test_idt_file_structure(self, idt_test_file):
        """Test that IDT file has the expected structure."""

        df = pd.read_excel(idt_test_file)

        # Basic structure checks
        assert len(df) > 0, "IDT file should have sequences"
        assert "Pool name" in df.columns, "Should have Pool name column"
        assert "Sequence" in df.columns, "Should have Sequence column"

        # Check for expected prefixes
        pool_names = df["Pool name"].dropna()
        fwd_count = sum(
            name.startswith("ExfoProbe_hIBD_noLink_P1") for name in pool_names
        )
        rev_count = sum(
            name.startswith("ExfoProbe_hIBD_noLink_P2") for name in pool_names
        )

        assert fwd_count > 0, "Should have forward primers"
        assert rev_count > 0, "Should have reverse primers"

        # Check sequences are not empty
        sequences = df["Sequence"].dropna()
        assert len(sequences) > 0, "Should have primer sequences"

        # Check sequence format (should contain DNA bases)
        sample_seq = sequences.iloc[0]
        dna_bases = set("ATCGN")
        seq_bases = set(sample_seq.upper())
        assert len(seq_bases & dna_bases) > 0, "Sequences should contain DNA bases"

        print(
            f"✓ IDT file validated: {len(df)} total sequences, {fwd_count} forward, {rev_count} reverse"
        )

    def test_recover_cli_missing_args(self):
        """Test recover CLI with missing required arguments."""
        import subprocess
        import sys

        # Test CLI with missing arguments
        result = subprocess.run(
            [sys.executable, "-m", "folitools.primer_selection.cli", "recover"],
            capture_output=True,
            text=True,
            cwd="/home/ubuntu/dev/folitools",
        )

        # Should exit with error code due to missing required arguments
        assert result.returncode != 0, "Missing args should cause non-zero exit"

        # Should show error message (might be in stdout or stderr)
        error_text = result.stdout + result.stderr
        assert len(error_text) > 0, "Should show error message for missing args"
        assert "order-excel" in error_text.lower(), (
            "Error should mention missing order-excel parameter"
        )
