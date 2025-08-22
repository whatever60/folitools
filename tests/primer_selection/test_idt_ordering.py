#!/usr/bin/env python3
"""
Test module for IDT ordering functionality in _04_make_excel.py
"""

import pytest
import pandas as pd
import tempfile
import os

from folitools.primer_selection._04_make_excel import summary


class TestIDTOrdering:
    """Test class for IDT ordering functionality."""

    def create_test_data(self, num_genes=5):
        """Create test data files for testing."""
        temp_dir = tempfile.mkdtemp()

        gene_names = [f"GENE_{i:03d}" for i in range(num_genes)]
        groups = [f"group_{i % 3}" for i in range(num_genes)]

        # Gene info file
        gene_data = pd.DataFrame({"gene": gene_names, "group": groups})
        gene_file = os.path.join(temp_dir, "genes.tsv")
        gene_data.to_csv(gene_file, sep="\t", index=False)

        # Primer selection file
        selection_data = pd.DataFrame(
            {
                "transcriptID": gene_names,
                "primerIndex": [f"{gene}_0" for gene in gene_names],
            }
        )
        selection_file = os.path.join(temp_dir, "selection.tsv")
        selection_data.to_csv(selection_file, sep="\t", index=False)

        # Primer info file
        primer_info_data = pd.DataFrame(
            {
                "geneSymbol": gene_names,
                "geneID": [f"ENSG{i:011d}" for i in range(num_genes)],
                "amplicon_index": [f"{gene}_0" for gene in gene_names],
                "L_seq": [f"ATGCCCGAGTC{i:03d}" for i in range(num_genes)],
                "R_seq": [f"TCAGGAAGCT{i:03d}" for i in range(num_genes)],
                "pool": [f"pool{i % 2}" for i in range(num_genes)],
            }
        )
        primer_info_file = os.path.join(temp_dir, "primer_info.tsv")
        primer_info_data.to_csv(primer_info_file, sep="\t", index=False)

        return temp_dir, gene_file, selection_file, primer_info_file

    def test_basic_idt_ordering(self):
        """Test basic IDT ordering functionality with small number of primers."""
        temp_dir, gene_file, selection_file, primer_info_file = self.create_test_data(5)

        try:
            output_excel = os.path.join(temp_dir, "summary.xlsx")
            output_idt = os.path.join(temp_dir, "idt_order.xlsx")

            # Run summary function
            summary(
                input_=gene_file,
                primer_selection=selection_file,
                primer_info=primer_info_file,
                output=output_excel,
                has_linker=False,
                output_idt_order=output_idt,
                idt_pool_prefix="TEST_POOL",
            )

            # Check that IDT file was created
            assert os.path.exists(output_idt), "IDT ordering file should be created"

            # Read and validate IDT file
            idt_data = pd.read_excel(output_idt)

            # Should have 2 primers per gene (forward + reverse)
            expected_primers = 5 * 2
            assert len(idt_data) == expected_primers, (
                f"Expected {expected_primers} primers"
            )

            # Check required columns
            assert "Pool name" in idt_data.columns, "Should have 'Pool name' column"
            assert "Sequence" in idt_data.columns, "Should have 'Sequence' column"

            # Check pool name (single pool, no suffix)
            unique_pools = idt_data["Pool name"].unique()
            assert len(unique_pools) == 1, "Should have single pool for small dataset"
            assert unique_pools[0] == "TEST_POOL", "Pool name should match prefix"

            # Check that sequences are properly formatted
            sequences = idt_data["Sequence"].tolist()
            for seq in sequences:
                assert isinstance(seq, str), "Sequences should be strings"
                assert len(seq) > 20, "Sequences should include prefixes"

        finally:
            import shutil

            shutil.rmtree(temp_dir)

    def test_multiple_pools(self):
        """Test IDT ordering with multiple pools (>384 primers)."""
        # Create 200 genes = 400 primers (should create 2 pools)
        temp_dir, gene_file, selection_file, primer_info_file = self.create_test_data(
            200
        )

        try:
            output_excel = os.path.join(temp_dir, "summary.xlsx")
            output_idt = os.path.join(temp_dir, "idt_order.xlsx")

            # Run summary function
            summary(
                input_=gene_file,
                primer_selection=selection_file,
                primer_info=primer_info_file,
                output=output_excel,
                has_linker=False,
                output_idt_order=output_idt,
                idt_pool_prefix="MULTI_POOL",
            )

            # Read IDT file
            idt_data = pd.read_excel(output_idt)

            # Check total primers
            assert len(idt_data) == 400, "Should have 400 primers for 200 genes"

            # Check pool distribution
            pool_counts = idt_data["Pool name"].value_counts()
            assert len(pool_counts) == 2, "Should have 2 pools"
            assert pool_counts["MULTI_POOL-p1"] == 384, (
                "First pool should have 384 primers"
            )
            assert pool_counts["MULTI_POOL-p2"] == 16, (
                "Second pool should have 16 primers"
            )

        finally:
            import shutil

            shutil.rmtree(temp_dir)

    def test_exact_384_primers(self):
        """Test with exactly 384 primers (should not add suffix)."""
        # 192 genes = 384 primers
        temp_dir, gene_file, selection_file, primer_info_file = self.create_test_data(
            192
        )

        try:
            output_excel = os.path.join(temp_dir, "summary.xlsx")
            output_idt = os.path.join(temp_dir, "idt_order.xlsx")

            summary(
                input_=gene_file,
                primer_selection=selection_file,
                primer_info=primer_info_file,
                output=output_excel,
                has_linker=False,
                output_idt_order=output_idt,
                idt_pool_prefix="EXACT_POOL",
            )

            idt_data = pd.read_excel(output_idt)

            # Check that there's exactly one pool with no suffix
            unique_pools = idt_data["Pool name"].unique()
            assert len(unique_pools) == 1, "Should have single pool"
            assert unique_pools[0] == "EXACT_POOL", (
                "Should not have suffix for 384 primers"
            )
            assert len(idt_data) == 384, "Should have exactly 384 primers"

        finally:
            import shutil

            shutil.rmtree(temp_dir)

    def test_with_linker(self):
        """Test IDT ordering with linker sequences."""
        temp_dir, gene_file, selection_file, primer_info_file = self.create_test_data(3)

        try:
            output_excel = os.path.join(temp_dir, "summary.xlsx")
            output_idt = os.path.join(temp_dir, "idt_order.xlsx")

            summary(
                input_=gene_file,
                primer_selection=selection_file,
                primer_info=primer_info_file,
                output=output_excel,
                has_linker=True,  # Test with linker
                output_idt_order=output_idt,
                idt_pool_prefix="LINKER_POOL",
            )

            idt_data = pd.read_excel(output_idt)

            # Check that sequences are longer (due to linker)
            sequences = idt_data["Sequence"].tolist()
            for seq in sequences:
                # Linker sequences are longer than non-linker
                assert len(seq) > 30, "Sequences with linker should be longer"

        finally:
            import shutil

            shutil.rmtree(temp_dir)

    def test_no_idt_output(self):
        """Test that function works normally when no IDT output is requested."""
        temp_dir, gene_file, selection_file, primer_info_file = self.create_test_data(3)

        try:
            output_excel = os.path.join(temp_dir, "summary.xlsx")

            # Run without IDT output
            result_df = summary(
                input_=gene_file,
                primer_selection=selection_file,
                primer_info=primer_info_file,
                output=output_excel,
                has_linker=False,
                output_idt_order=None,  # No IDT output
                idt_pool_prefix="UNUSED",
            )

            # Should work normally and return DataFrame
            assert isinstance(result_df, pd.DataFrame), "Should return DataFrame"
            assert len(result_df) == 3, "Should have 3 genes"

        finally:
            import shutil

            shutil.rmtree(temp_dir)

    def test_gene_pair_consistency(self):
        """Test that forward and reverse primers for same gene are in same pool."""
        temp_dir, gene_file, selection_file, primer_info_file = self.create_test_data(
            200
        )

        try:
            output_excel = os.path.join(temp_dir, "summary.xlsx")
            output_idt = os.path.join(temp_dir, "idt_order.xlsx")

            summary(
                input_=gene_file,
                primer_selection=selection_file,
                primer_info=primer_info_file,
                output=output_excel,
                has_linker=False,
                output_idt_order=output_idt,
                idt_pool_prefix="CONSISTENCY_TEST",
            )

            idt_data = pd.read_excel(output_idt)

            # Check that pairs are in same pool
            for i in range(0, len(idt_data), 2):
                if i + 1 < len(idt_data):
                    pool1 = idt_data.iloc[i]["Pool name"]
                    pool2 = idt_data.iloc[i + 1]["Pool name"]
                    assert pool1 == pool2, (
                        f"Primer pair at positions {i}, {i + 1} should be in same pool"
                    )

        finally:
            import shutil

            shutil.rmtree(temp_dir)


if __name__ == "__main__":
    pytest.main([__file__])
