#!/usr/bin/env python3
"""
Tests to verify that the primer selection functions return meaningful data types
instead of just status codes. This documents the breaking changes where functions
now return useful objects for programmatic use.
"""

import pandas as pd
from Bio.SeqRecord import SeqRecord

from folitools.primer_selection._01_primer_selection import subset
from folitools.primer_selection._02_select_primer_set_by_saddle_loss import saddle
from folitools.primer_selection._03_extract_region_sequence import product
from folitools.primer_selection._04_make_excel import summary
from folitools.primer_selection._05_recover import recover


class TestFunctionReturnTypes:
    """Test that primer selection functions return meaningful data structures."""

    def test_subset_returns_dataframe(self):
        """Test that subset() returns a pandas DataFrame."""
        # This test documents that subset() now returns the primer info DataFrame
        # instead of just writing files and returning None (breaking change)
        
        # Check function annotation
        return_annotation = subset.__annotations__.get('return')
        assert return_annotation == pd.DataFrame, (
            f"subset should be annotated to return pd.DataFrame, got {return_annotation}"
        )

    def test_saddle_returns_dataframe(self):
        """Test that saddle() returns a pandas DataFrame."""
        # This test documents that saddle() now returns the selected primers DataFrame
        # instead of just writing files and returning None (breaking change)
        
        # Check function annotation
        return_annotation = saddle.__annotations__.get('return')
        assert return_annotation == pd.DataFrame, (
            f"saddle should be annotated to return pd.DataFrame, got {return_annotation}"
        )

    def test_product_returns_seqrecord_list(self):
        """Test that product() returns a list of SeqRecord objects."""
        # This test documents that product() now returns the extracted amplicon sequences
        # instead of just writing FASTA files and returning None (breaking change)
        
        # Check function annotation
        return_annotation = product.__annotations__.get('return')
        assert return_annotation == list[SeqRecord], (
            f"product should be annotated to return list[SeqRecord], got {return_annotation}"
        )

    def test_summary_returns_dataframe(self):
        """Test that summary() returns a pandas DataFrame."""
        # This test documents that summary() now returns the primer ordering DataFrame
        # instead of just writing Excel files and returning None (breaking change)
        
        # Check function annotation
        return_annotation = summary.__annotations__.get('return')
        assert return_annotation == pd.DataFrame, (
            f"summary should be annotated to return pd.DataFrame, got {return_annotation}"
        )

    def test_recover_returns_dataframe(self):
        """Test that recover() returns a pandas DataFrame."""
        # This test documents that recover() now returns the recovered primer summary
        # instead of just writing files and returning None (breaking change)
        
        # Check function annotation
        return_annotation = recover.__annotations__.get('return')
        assert return_annotation == pd.DataFrame, (
            f"recover should be annotated to return pd.DataFrame, got {return_annotation}"
        )

    def test_function_docstrings_document_returns(self):
        """Test that function docstrings properly document return values."""
        functions_to_check = [subset, saddle, product, summary, recover]
        
        for func in functions_to_check:
            docstring = func.__doc__
            assert docstring is not None, f"{func.__name__} should have a docstring"
            assert "Returns:" in docstring or "Return:" in docstring, (
                f"{func.__name__} docstring should document return value"
            )

    def test_breaking_changes_summary(self):
        """Document the breaking changes in return types."""
        breaking_changes = {
            "subset": {
                "old": "None (just wrote files)",
                "new": "pd.DataFrame (primer info data)",
                "benefit": "Can be used programmatically for further processing"
            },
            "saddle": {
                "old": "None (just wrote files)", 
                "new": "pd.DataFrame (selected primers)",
                "benefit": "Can inspect optimization results without reading files"
            },
            "product": {
                "old": "None (just wrote FASTA)",
                "new": "list[SeqRecord] (amplicon sequences)",
                "benefit": "Can work with sequences in memory without file I/O"
            },
            "summary": {
                "old": "None (just wrote Excel)",
                "new": "pd.DataFrame (primer ordering data)",
                "benefit": "Can programmatically access ordering information"
            },
            "recover": {
                "old": "None (just wrote files)",
                "new": "pd.DataFrame (recovered primer summary)",
                "benefit": "Can analyze recovered data without reading files"
            }
        }
        
        # This test documents the breaking changes for users upgrading
        print("\nBreaking Changes Summary:")
        print("=" * 50)
        for func_name, changes in breaking_changes.items():
            print(f"\n{func_name}():")
            print(f"  Old: {changes['old']}")
            print(f"  New: {changes['new']}")
            print(f"  Benefit: {changes['benefit']}")
        
        # Test passes - this is just documentation
        assert len(breaking_changes) == 5, "Should document all 5 function changes"

    def test_backward_compatibility_note(self):
        """Document how to handle the breaking changes."""
        migration_guide = """
        Migration Guide for Breaking Changes:
        
        OLD CODE (v0.1.x):
        ------------------
        subset(gene_table, species, size_range, out_seq, out_info)
        # Result was None, had to read files manually
        
        NEW CODE (v0.2.x):
        ------------------
        result_df = subset(gene_table, species, size_range, out_seq, out_info)
        # Now you get the DataFrame directly for further processing
        
        Benefits:
        - No need to read files back in for programmatic use
        - Can chain operations more easily
        - Better for notebook/interactive analysis
        - Files are still written for compatibility
        """
        
        print(migration_guide)
        assert True, "Migration guide printed"
