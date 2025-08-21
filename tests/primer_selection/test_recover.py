"""Tests for the _05_recover module."""

import logging

import pandas as pd
import pytest

from folitools.primer_selection._05_recover import (
    _read_idt_excel,
    _classify_primers,
    _sanity_check_primer_duplicates,
    _sanity_check_sequence_uniqueness,
    _sanity_check_seqkit_patterns,
    _sanity_check_multi_location_binding,
    _sanity_check_multi_gene_binding,
    _sanity_check_sequence_lengths,
    _sanity_check_genes_per_primer_pair,
    _sanity_check_amplicons_per_pair,
    _sanity_check_cross_pool_amplicons,
    _sanity_check_primer_pair_relationships,
    _enrich_locate_df,
    _build_amplicons,
    _group_pairs,
    _build_summary_df,
)
from folitools.primer_selection.utils import get_prefixes


class TestReadIdtExcel:
    """Test the _read_idt_excel function."""

    def test_read_idt_excel_basic(self, sample_idt_excel):
        """Test basic IDT Excel reading functionality."""
        df = _read_idt_excel(sample_idt_excel)

        assert len(df) == 4
        assert list(df.columns) == ["pool", "sequence"]
        assert all(df["sequence"].str.contains(r"^[ACGTN]+$"))
        assert all(df["sequence"].str.isupper())

    def test_read_idt_excel_robust_false(self, sample_idt_excel):
        """Test non-robust reading mode."""
        df = _read_idt_excel(sample_idt_excel, robust=False)

        assert len(df) == 4
        assert list(df.columns) == ["pool", "sequence"]


class TestClassifyPrimers:
    """Test the _classify_primers function."""

    def test_classify_primers_success(self, sample_idt_excel):
        """Test successful primer classification."""
        df_idt = _read_idt_excel(sample_idt_excel)
        # Use the correct prefixes for has_linker=True
        fwd_prefix, rev_prefix = get_prefixes(has_linker=True)

        result = _classify_primers(df_idt, fwd_prefix, rev_prefix)

        assert len(result) == 4
        expected_cols = ["pool", "primer_type", "primer_seq_full", "primer_seq"]
        assert list(result.columns) == expected_cols
        assert set(result["primer_type"]) == {"fwd", "rev"}

    def test_classify_primers_unmatched(self, sample_idt_excel):
        """Test primer classification with unmatched sequences."""
        df_idt = _read_idt_excel(sample_idt_excel)
        # Use wrong prefixes to trigger ValueError
        fwd_prefix = "WRONG_PREFIX"
        rev_prefix = "WRONG_PREFIX"

        with pytest.raises(ValueError, match="Some primers do not match"):
            _classify_primers(df_idt, fwd_prefix, rev_prefix)


class TestSanityChecks:
    """Test all sanity check functions."""

    def test_sanity_check_primer_duplicates_no_duplicates(
        self, sample_primer_info, caplog
    ):
        """Test sanity check with no duplicate primers."""
        with caplog.at_level(logging.WARNING):
            _sanity_check_primer_duplicates(sample_primer_info)

        # Should not log any warnings for unique primers
        assert len(caplog.records) == 0

    def test_sanity_check_primer_duplicates_with_duplicates(self, caplog):
        """Test sanity check with duplicate primers."""
        # Create DataFrame with duplicates
        primer_info_dup = pd.DataFrame(
            {
                "primer_seq": ["ATGAAGACGGCATC", "ATGAAGGCTTCC", "ATGAAGACGGCATC"],
                "pool": ["Pool1", "Pool2", "Pool3"],
            }
        )

        with caplog.at_level(logging.WARNING):
            _sanity_check_primer_duplicates(primer_info_dup)

        assert len(caplog.records) > 0
        assert "Duplicate primer gene-specific sequences detected" in caplog.text
        assert "ATGAAGACGGCATC" in caplog.text

    def test_sanity_check_sequence_uniqueness_unique(self, caplog):
        """Test sequence uniqueness check with unique sequences."""
        primer_seqs = ["ATGAAGACGGCATC", "ATGAAGGCTTCC", "ATCCGATGCATTAG"]

        with caplog.at_level(logging.INFO):
            _sanity_check_sequence_uniqueness(primer_seqs)

        assert "✅ All primer sequences are unique" in caplog.text

    def test_sanity_check_sequence_uniqueness_duplicates(self, caplog):
        """Test sequence uniqueness check with duplicates."""
        primer_seqs = ["ATGAAGACGGCATC", "ATGAAGGCTTCC", "ATGAAGACGGCATC"]

        with caplog.at_level(logging.INFO):
            _sanity_check_sequence_uniqueness(primer_seqs)

        assert "Duplicate primer sequences found" in caplog.text
        assert "ATGAAGACGGCATC: 2 occurrences" in caplog.text

    def test_sanity_check_seqkit_patterns_unique(self, sample_locate_df, caplog):
        """Test seqkit patterns check with unique patterns."""
        with caplog.at_level(logging.INFO):
            _sanity_check_seqkit_patterns(sample_locate_df)

        assert "✅ All seqkit locate patterns are unique" in caplog.text

    def test_sanity_check_multi_location_binding_none(self, caplog):
        """Test multi-location binding check with no multi-location primers."""
        locate_df = pd.DataFrame(
            {
                "primer_seq": ["ATGAAGACGGCATC", "ATGAAGGCTTCC"],
                "transcript_id": ["ENSMUST00000000001.4", "ENSMUST00000000002.4"],
            }
        )

        with caplog.at_level(logging.INFO):
            _sanity_check_multi_location_binding(locate_df)

        assert (
            "✅ No primers bind to multiple locations on the same transcript"
            in caplog.text
        )

    def test_sanity_check_multi_gene_binding_none(self, caplog):
        """Test multi-gene binding check with no multi-gene primers."""
        locate_df = pd.DataFrame(
            {
                "primer_seq": ["ATGAAGACGGCATC", "ATGAAGGCTTCC"],
                "gene_id": ["ENSMUSG00000000001.4", "ENSMUSG00000000002.4"],
            }
        )

        with caplog.at_level(logging.INFO):
            _sanity_check_multi_gene_binding(locate_df)

        assert "✅ All primers bind to only one gene each" in caplog.text

    def test_sanity_check_sequence_lengths_match(self, caplog):
        """Test sequence length check with matching lengths."""
        locate_df = pd.DataFrame(
            {
                "primer_seq": ["ATGAAGACGGCATC", "ATGAAGGCTTCC"],
                "start": [100, 200],
                "end": [113, 211],  # 14 and 12 bases respectively
            }
        )

        with caplog.at_level(logging.INFO):
            _sanity_check_sequence_lengths(locate_df)

        assert (
            "✅ All primer sequence lengths match their locate positions" in caplog.text
        )

    def test_sanity_check_genes_per_primer_pair(self, sample_amplicon_df, caplog):
        """Test genes per primer pair analysis."""
        with caplog.at_level(logging.INFO):
            _sanity_check_genes_per_primer_pair(sample_amplicon_df)

        assert "Analyzing genes amplified per primer pair" in caplog.text
        assert "1 gene(s): 2 primer pairs" in caplog.text

    def test_sanity_check_amplicons_per_pair(self, sample_amplicon_df, caplog):
        """Test amplicons per pair analysis."""
        with caplog.at_level(logging.INFO):
            _sanity_check_amplicons_per_pair(sample_amplicon_df)

        assert "Analyzing amplicons formed per primer pair" in caplog.text
        assert "Summary: min=1, max=1 amplicons per pair" in caplog.text

    def test_sanity_check_cross_pool_amplicons(self, sample_amplicon_df, caplog):
        """Test cross-pool amplicon analysis."""
        with caplog.at_level(logging.INFO):
            _sanity_check_cross_pool_amplicons(sample_amplicon_df)

        assert "Checking for cross-pool amplicons" in caplog.text
        assert "Total amplicons: 2" in caplog.text
        assert "Cross pool (different pools): 2 (100.0%)" in caplog.text

    def test_sanity_check_primer_pair_relationships_good(self, caplog):
        """Test primer pair relationships with good one-to-one mapping."""
        grouped_df = pd.DataFrame(
            {"primer_seq_fwd": ["SEQ1", "SEQ2"], "primer_seq_rev": ["SEQ3", "SEQ4"]}
        )

        with caplog.at_level(logging.INFO):
            _sanity_check_primer_pair_relationships(grouped_df)

        assert "✅ One-to-one primer relationships maintained" in caplog.text


class TestHelperFunctions:
    """Test helper functions."""

    def test_enrich_locate_df(self, sample_locate_df, sample_primer_info):
        """Test enriching locate DataFrame with primer info."""
        result = _enrich_locate_df(sample_locate_df, sample_primer_info)

        expected_cols = [
            "primer_seq",
            "primer_seq_full",
            "primer_type",
            "transcript_id",
            "gene_id",
            "gene_symbol",
            "start",
            "end",
            "strand",
            "pool",
        ]
        assert list(result.columns) == expected_cols
        assert len(result) == 2
        assert "Gnai3" in result["gene_symbol"].values

    def test_build_amplicons_empty_result(self):
        """Test building amplicons with input that produces no amplicons."""
        # Create locate_final with only forward primers (no pairs possible)
        locate_final = pd.DataFrame(
            {
                "transcript_id": ["ENSMUST00000000001.4"],
                "primer_type": ["fwd"],
                "primer_seq": ["ATGAAGACGGCATC"],
                "gene_id": ["ENSMUSG00000000001.4"],
                "gene_symbol": ["Gnai3"],
                "start": [100],
                "end": [113],
                "strand": ["+"],
                "pool": ["Pool1"],
            }
        )

        with pytest.raises(ValueError, match="No amplicons were formed"):
            _build_amplicons(locate_final)

    def test_group_pairs_empty_input(self):
        """Test grouping pairs with empty input."""
        empty_df = pd.DataFrame()
        result = _group_pairs(empty_df)

        expected_cols = [
            "primer_seq_fwd",
            "primer_seq_rev",
            "transcript_id",
            "gene_id",
            "gene_symbol",
            "start_up",
            "end_up",
            "start_down",
            "end_down",
            "pool_fwd",
            "pool_rev",
            "num_transcripts",
            "num_genes",
            "transcript_id_all",
        ]
        assert list(result.columns) == expected_cols
        assert len(result) == 0

    def test_build_summary_df_empty_input(self):
        """Test building summary DataFrame with empty input."""
        empty_df = pd.DataFrame()
        fwd_prefix = "CCTACACGACGCTCTTCCGATCTNNNNNNACATCA"
        rev_prefix = "TCAGACGTGTGCTCTTCCGATCTNNNNNNATAGTTT"

        result = _build_summary_df(empty_df, fwd_prefix, rev_prefix)

        expected_cols = [
            "Group",
            "geneSymbol",
            "geneID",
            "Chosen Index",
            "amplicon_index",
            "L_seq",
            "R_seq",
            "primer_sequence_to_order_forward",
            "primer_sequence_to_order_reverse",
            "L_pool",
            "R_pool",
            "multi_mapping",
        ]
        assert list(result.columns) == expected_cols
        assert len(result) == 0
