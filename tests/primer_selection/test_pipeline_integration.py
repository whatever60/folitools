#!/usr/bin/env python3
"""
Integration tests for primer selection pipeline.
"""

import pytest
import pandas as pd
import tempfile
import shutil
from pathlib import Path

from folitools.primer_selection._01_primer_selection import subset
from folitools.primer_selection._02_select_primer_set_by_saddle_loss import saddle
from folitools.primer_selection.cli import workflow


class TestPipelineIntegration:
    """Integration tests for primer selection pipeline."""

    @pytest.fixture
    def mouse_genes_file(self):
        """Use the test mouse genes file."""
        genes_file = Path(__file__).parent / "test_data" / "mouse_genes.tsv"
        if not genes_file.exists():
            pytest.skip(f"Mouse genes file not found: {genes_file}")
        return genes_file

    @pytest.fixture
    def temp_output_dir(self):
        """Create temporary output directory."""
        temp_dir = Path(tempfile.mkdtemp(prefix="test_pipeline_"))
        yield temp_dir
        # Cleanup
        if temp_dir.exists():
            shutil.rmtree(temp_dir)

    def test_step1_subset(self, mouse_genes_file, temp_output_dir):
        """Test subset step."""
        output_primer_sequence = temp_output_dir / "primer_sequence.tsv"
        output_primer_info = temp_output_dir / "primer_info.tsv"

        result_df = subset(
            gene_table_file=mouse_genes_file,
            species="mouse",
            amplicon_size_range=(320, 380),
            output_primer_sequence=output_primer_sequence,
            output_primer_info=output_primer_info,
        )

        assert result_df is not None, "Subset should return a DataFrame"
        assert not result_df.empty, "Subset should return non-empty DataFrame"
        assert output_primer_sequence.exists(), "Primer sequence file should be created"
        assert output_primer_info.exists(), "Primer info file should be created"

        # Check file contents
        primer_seq_df = pd.read_csv(output_primer_sequence, sep="\t")
        primer_info_df = pd.read_csv(output_primer_info, sep="\t")

        assert len(primer_seq_df) > 0, "Should find primer sequences"
        assert len(primer_info_df) > 0, "Should find primer info"

        print(
            f"✓ Subset: Found {len(primer_seq_df)} primer sequences for {len(primer_info_df)} amplicons"
        )

    def test_step2_saddle(self, mouse_genes_file, temp_output_dir):
        """Test saddle step."""
        # First run subset to get input files
        primer_sequence_file = temp_output_dir / "primer_sequence.tsv"
        primer_info_file = temp_output_dir / "primer_info.tsv"

        subset_result = subset(
            gene_table_file=mouse_genes_file,
            species="mouse",
            amplicon_size_range=(320, 380),
            output_primer_sequence=primer_sequence_file,
            output_primer_info=primer_info_file,
        )

        assert subset_result is not None, "Subset should return a DataFrame"
        assert not subset_result.empty, "Subset should return non-empty DataFrame"

        # Run saddle
        output_selected = temp_output_dir / "selected_primers.tsv"
        output_loss = temp_output_dir / "saddle_loss.txt"

        result_df = saddle(
            input_=primer_sequence_file,
            output=output_selected,
            output_loss=output_loss,
            num_cycles_anneal=5,  # Reduced for testing
            random_seed=42,
            background_fasta=None,
        )

        assert result_df is not None, "Saddle should return a DataFrame"
        assert not result_df.empty, "Saddle should return non-empty DataFrame"
        assert output_selected.exists(), "Selected primers file should be created"
        assert output_loss.exists(), "Loss file should be created"

        # Check results
        selected_df = pd.read_csv(output_selected, sep="\t")
        assert len(selected_df) > 0, "Should select some primers"

        print(f"✓ Saddle: Selected {len(selected_df)} optimal primer pairs")

    def test_workflow_command(self, mouse_genes_file, temp_output_dir):
        """Test the complete workflow command."""

        # Run workflow with subset of genes for faster testing
        small_genes_file = temp_output_dir / "small_genes.tsv"
        genes_df = pd.read_csv(mouse_genes_file, sep="\t")
        small_genes_df = genes_df.head(20)  # Use only 20 genes for speed
        small_genes_df.to_csv(small_genes_file, sep="\t", index=False)

        exit_code = workflow(
            input_=small_genes_file,
            species="mouse",
            amplicon_size_range=(320, 380),
            output_dir=temp_output_dir,
            num_cycles_anneal=3,  # Very reduced for testing
        )

        assert exit_code == 0, "Workflow should complete successfully"

        # Check all expected output files
        expected_files = [
            "primer_sequence.tsv",
            "primer_info.tsv",
            "primer_sequence_selected.tsv",
            "primer_sequence_selected_loss.txt",
            "amplicons.fasta",
            "primer_to_order.xlsx",
        ]

        for filename in expected_files:
            filepath = temp_output_dir / filename
            assert filepath.exists(), f"Expected output file {filename} should exist"

        # Check results
        selected_df = pd.read_csv(
            temp_output_dir / "primer_sequence_selected.tsv", sep="\t"
        )
        assert len(selected_df) > 0, "Workflow should select primers"

        print(
            f"✓ Workflow: Generated {len(expected_files)} output files with {len(selected_df)} selected primers"
        )

    def test_gene_file_format(self, mouse_genes_file):
        """Test gene file format."""
        assert mouse_genes_file.exists(), "Gene file should exist"

        df = pd.read_csv(mouse_genes_file, sep="\t")
        assert len(df) > 0, "Gene file should not be empty"

        # Check required columns
        required_cols = ["gene"]  # Actual column name in the gene file
        for col in required_cols:
            assert col in df.columns, f"Gene file should have {col} column"

        # Check for group information if present
        if "group" in df.columns:
            groups = df["group"].value_counts()
            print(f"✓ Gene file: {len(df)} genes across {len(groups)} groups")
        else:
            print(f"✓ Gene file: {len(df)} genes loaded successfully")
