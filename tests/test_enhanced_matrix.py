"""Test the enhanced read_counts function with gene co-occurrence clustering."""

import pandas as pd
import polars as pl
import tempfile
import os

from folitools.get_matrix import (
    consolidate_read_id, 
    per_read_to_edges, 
    edges_to_graph, 
    extract_all_genes,
    process_count_file,
    read_counts
)


def test_gene_cooccurrence_clustering():
    """Test that genes co-occurring in reads are properly clustered."""
    # Create test data where genes co-occur in reads
    data = {
        "read_id": ["read1", "read1", "read2", "read2", "read3", "read3", "read4"],
        "contig": ["chr1", "chr1", "chr2", "chr2", "chr3", "chr3", "chr4"],
        "position": [100, 200, 300, 400, 500, 600, 700],
        "gene": ["GENE_A", "GENE_B", "GENE_A", "GENE_B", "GENE_C", "GENE_D", "GENE_E"],
        "umi": ["AAA", "BBB", "CCC", "DDD", "EEE", "FFF", "GGG"],
        "umi_count": [1, 1, 1, 1, 1, 1, 1],
        "final_umi": ["AAA", "BBB", "CCC", "DDD", "EEE", "FFF", "GGG"],
        "final_umi_count": [1, 1, 1, 1, 1, 1, 1],
        "unique_id": ["1", "2", "3", "4", "5", "6", "7"]
    }
    
    df_test = pd.DataFrame(data)
    
    # Create temporary test file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
        df_test.to_csv(f.name, sep='\t', index=False)
        temp_file = f.name
    
    try:
        # Process the file
        result = process_count_file(temp_file)
        
        # Verify clustering: GENE_A and GENE_B should be clustered together
        # since they co-occur in read1 and read2
        gene_names = list(result.index)
        
        # There should be a cluster containing both GENE_A and GENE_B
        found_ab_cluster = any("GENE_A" in gene and "GENE_B" in gene for gene in gene_names)
        assert found_ab_cluster, f"GENE_A and GENE_B should be clustered together. Found genes: {gene_names}"
        
        # GENE_C and GENE_D should be clustered together (co-occur in read3)
        found_cd_cluster = any("GENE_C" in gene and "GENE_D" in gene for gene in gene_names)
        assert found_cd_cluster, f"GENE_C and GENE_D should be clustered together. Found genes: {gene_names}"
        
        # GENE_E should be isolated (only appears in read4)
        found_e_isolated = "GENE_E" in gene_names
        assert found_e_isolated, f"GENE_E should appear as isolated gene. Found genes: {gene_names}"
        
        # Total UMI count should be preserved
        total_umis = result.sum()
        assert total_umis == 7, f"Expected 7 total UMIs, got {total_umis}"
        
        print("✓ Gene co-occurrence clustering test passed")
        
    finally:
        # Clean up
        os.unlink(temp_file)


def test_graph_building_functions():
    """Test the individual graph building functions."""
    # Create test dataframe
    data = {
        "read_id": ["read1", "read1", "read2", "read2", "read3"],
        "gene": ["GENE_A", "GENE_B", "GENE_A", "GENE_C", "GENE_D"]
    }
    df = pl.DataFrame(data).lazy()  # Convert to LazyFrame
    
    # Test per-read building
    df_per_read = consolidate_read_id(df)
    df_per_read_collected = df_per_read.collect()
    assert df_per_read_collected.shape[0] == 2, "Should have 2 reads with multiple genes"
    
    # Test edge generation
    edges_df = per_read_to_edges(df_per_read)
    edges_df_collected = edges_df.collect()
    assert edges_df_collected.shape[0] >= 1, "Should have at least one edge"
    
    # Test gene extraction
    all_genes = extract_all_genes(df.collect())  # extract_all_genes still takes DataFrame
    assert len(all_genes) == 4, f"Should have 4 unique genes, got {len(all_genes)}"
    assert set(all_genes) == {"GENE_A", "GENE_B", "GENE_C", "GENE_D"}
    
    # Test graph creation
    graph = edges_to_graph(edges_df, all_genes=all_genes)
    assert graph.number_of_nodes() == 4, "Graph should have 4 nodes"
    
    print("✓ Graph building functions test passed")


def test_read_counts_integration():
    """Test the full read_counts function with clustering."""
    # Create multiple test files
    test_files = []
    
    for sample_idx in range(2):
        data = {
            "read_id": [f"read{sample_idx}_1", f"read{sample_idx}_1", f"read{sample_idx}_2", f"read{sample_idx}_2", f"read{sample_idx}_3"],
            "contig": ["chr1"] * 5,
            "position": [100 + i * 100 for i in range(5)],
            "gene": ["GENE_A", "GENE_B", "GENE_A", "GENE_B", "GENE_C"],  # GENE_A and GENE_B co-occur
            "umi": [f"UMI_{i}" for i in range(5)],
            "umi_count": [1] * 5,
            "final_umi": [f"UMI_{i}" for i in range(5)],
            "final_umi_count": [1] * 5,
            "unique_id": [str(i) for i in range(5)]
        }
        
        df_test = pd.DataFrame(data)
        
        # Create temporary test file
        with tempfile.NamedTemporaryFile(mode='w', suffix=f'.sample{sample_idx}.tsv', delete=False) as f:
            df_test.to_csv(f.name, sep='\t', index=False)
            test_files.append(f.name)
    
    try:
        # Test read_counts function
        matrix = read_counts(test_files)
        
        # Check matrix dimensions
        assert matrix.shape[0] == 2, f"Should have 2 samples, got {matrix.shape[0]}"
        
        # Check that clustered genes appear in column names
        col_names = list(matrix.columns)
        has_cluster = any("|" in col for col in col_names)
        assert has_cluster, f"Should have clustered gene names with '|', got columns: {col_names}"
        
        print("✓ read_counts integration test passed")
        
    finally:
        # Clean up
        for f in test_files:
            os.unlink(f)


if __name__ == "__main__":
    test_graph_building_functions()
    test_gene_cooccurrence_clustering()
    test_read_counts_integration()
    print("\n✓ All enhanced read_counts tests passed!")
