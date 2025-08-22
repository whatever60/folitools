#!/usr/bin/env python3
"""
Test script to verify IDT ordering functionality.
"""

import pandas as pd
import tempfile
import os

# Add the source directory to Python path
import sys

sys.path.insert(0, "/home/ubuntu/dev/folitools/src")

from folitools.primer_selection._04_make_excel import summary


def create_test_data():
    """Create test data files for testing IDT ordering."""

    # Create temporary directory
    temp_dir = tempfile.mkdtemp()

    # Gene info file (input_)
    gene_data = pd.DataFrame(
        {
            "gene": ["ACTB", "GAPDH", "TNF", "IL1B", "TP53"],
            "group": [
                "housekeeping",
                "housekeeping",
                "cytokine",
                "cytokine",
                "tumor_suppressor",
            ],
        }
    )
    gene_file = os.path.join(temp_dir, "genes.tsv")
    gene_data.to_csv(gene_file, sep="\t", index=False)

    # Primer selection file
    selection_data = pd.DataFrame(
        {
            "transcriptID": ["ACTB", "GAPDH", "TNF", "IL1B", "TP53"],
            "primerIndex": ["ACTB_0", "GAPDH_0", "TNF_1", "IL1B_0", "TP53_2"],
        }
    )
    selection_file = os.path.join(temp_dir, "selection.tsv")
    selection_data.to_csv(selection_file, sep="\t", index=False)

    # Primer info file
    primer_info_data = pd.DataFrame(
        {
            "geneSymbol": ["ACTB", "GAPDH", "TNF", "IL1B", "TP53"],
            "geneID": [
                "ENSG00000075624",
                "ENSG00000111640",
                "ENSG00000232810",
                "ENSG00000125538",
                "ENSG00000141510",
            ],
            "amplicon_index": ["ACTB_0", "GAPDH_0", "TNF_1", "IL1B_0", "TP53_2"],
            "L_seq": [
                "ATGCCCGAGTC",
                "GGTGAAGGTCG",
                "CAGAGGGAAGC",
                "TTGTGGCTTCG",
                "CCACCAGCAGC",
            ],
            "R_seq": [
                "TCAGGAAGCT",
                "GCAGTGGTTCC",
                "CTGGTTATCGC",
                "AAGCTGATGG",
                "GGATGGTGGT",
            ],
            "pool": ["pool1", "pool1", "pool2", "pool2", "pool3"],
        }
    )
    primer_info_file = os.path.join(temp_dir, "primer_info.tsv")
    primer_info_data.to_csv(primer_info_file, sep="\t", index=False)

    return temp_dir, gene_file, selection_file, primer_info_file


def test_idt_ordering():
    """Test the IDT ordering functionality."""

    print("Creating test data...")
    temp_dir, gene_file, selection_file, primer_info_file = create_test_data()

    # Output files
    output_excel = os.path.join(temp_dir, "primer_summary.xlsx")
    output_idt = os.path.join(temp_dir, "idt_order.xlsx")

    print("Running summary function with IDT ordering...")
    try:
        summary(
            input_=gene_file,
            primer_selection=selection_file,
            primer_info=primer_info_file,
            output=output_excel,
            has_linker=False,
            output_idt_order=output_idt,
            idt_pool_prefix="TEST_POOL",
        )

        print("✓ Summary function completed successfully")

        # Check if IDT file was created
        if os.path.exists(output_idt):
            print("✓ IDT ordering file created")

            # Read and display the IDT file
            idt_data = pd.read_excel(output_idt)
            print(f"✓ IDT file contains {len(idt_data)} primers")
            print("\nIDT ordering preview:")
            print(idt_data.head(10))

            # Check pool names
            unique_pools = idt_data["Pool name"].unique()
            print(f"\nUnique pools: {list(unique_pools)}")

        else:
            print("✗ IDT ordering file was not created")

    except Exception as e:
        print(f"✗ Error running summary function: {e}")
        import traceback

        traceback.print_exc()

    finally:
        # Cleanup
        import shutil

        shutil.rmtree(temp_dir)
        print(f"\nCleaned up temporary directory: {temp_dir}")


if __name__ == "__main__":
    test_idt_ordering()
