#!/usr/bin/env python3
"""
Test script to verify IDT ordering functionality with exactly 384 primers.
"""

import pandas as pd
import tempfile
import os

# Add the source directory to Python path
import sys

sys.path.insert(0, "/home/ubuntu/dev/folitools/src")

from folitools.primer_selection._04_make_excel import summary


def test_exact_384_primers():
    """Test the IDT ordering functionality with exactly 384 primers (192 genes)."""

    print("Creating test data for exactly 384 primers...")
    temp_dir = tempfile.mkdtemp()

    # Create 192 genes (384 primers total) to test single pool behavior
    num_genes = 192
    gene_names = [f"GENE_{i:03d}" for i in range(num_genes)]
    groups = [f"group_{i % 3}" for i in range(num_genes)]  # 3 different groups

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
            "pool": [f"pool{i % 2}" for i in range(num_genes)],  # 2 different pools
        }
    )
    primer_info_file = os.path.join(temp_dir, "primer_info.tsv")
    primer_info_data.to_csv(primer_info_file, sep="\t", index=False)

    # Output files
    output_excel = os.path.join(temp_dir, "primer_summary.xlsx")
    output_idt = os.path.join(temp_dir, "idt_order.xlsx")

    print(f"Testing with {num_genes} genes ({num_genes * 2} primers)...")

    try:
        summary(
            input_=gene_file,
            primer_selection=selection_file,
            primer_info=primer_info_file,
            output=output_excel,
            has_linker=False,
            output_idt_order=output_idt,
            idt_pool_prefix="SINGLE_POOL_TEST",
        )

        print("✓ Summary function completed successfully")

        # Check if IDT file was created
        if os.path.exists(output_idt):
            print("✓ IDT ordering file created")

            # Read and analyze the IDT file
            idt_data = pd.read_excel(output_idt)
            print(f"✓ IDT file contains {len(idt_data)} primers")

            # Check pool names
            unique_pools = idt_data["Pool name"].unique()
            print(f"Unique pools: {list(unique_pools)}")

            # Check that there's no suffix when exactly 384 primers
            if len(unique_pools) == 1 and unique_pools[0] == "SINGLE_POOL_TEST":
                print("✓ Correct behavior: No suffix added for single pool")
            else:
                print("✗ Error: Unexpected pool naming for 384 primers")

            # Check pool distribution
            pool_counts = idt_data["Pool name"].value_counts().sort_index()
            print("Pool distribution:")
            for pool, count in pool_counts.items():
                print(f"  {pool}: {count} primers")

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
    test_exact_384_primers()
