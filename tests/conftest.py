"""Shared test fixtures for folitools tests."""

import tempfile
from pathlib import Path
from typing import Generator

import pandas as pd
import pytest


@pytest.fixture
def sample_idt_excel() -> Generator[Path, None, None]:
    """Create a temporary Excel file with sample IDT order data."""
    test_data = {
        'Pool': ['Pool1', 'Pool2', 'Pool3', 'Pool4'],
        'Sequence': [
            'CCTACACGACGCTCTTCCGATCTNNNNNNACATCAATGAAGACGGCATC',  # Forward primer
            'TCAGACGTGTGCTCTTCCGATCTNNNNNNATAGTTATGAAGGCTTCC',  # Reverse primer  
            'CCTACACGACGCTCTTCCGATCTNNNNNNACATCATCCGATGCATTAG',  # Different forward
            'TCAGACGTGTGCTCTTCCGATCTNNNNNNATAGTTATCGGAATGCTCC',  # Different reverse
        ]
    }
    
    df = pd.DataFrame(test_data)
    
    with tempfile.NamedTemporaryFile(suffix='.xlsx', delete=False) as tmp_file:
        df.to_excel(tmp_file.name, index=False)
        excel_path = Path(tmp_file.name)
    
    try:
        yield excel_path
    finally:
        excel_path.unlink(missing_ok=True)


@pytest.fixture
def sample_idt_excel_with_duplicates() -> Generator[Path, None, None]:
    """Create a temporary Excel file with duplicate primer sequences."""
    test_data = {
        'Pool': ['Pool1', 'Pool2', 'Pool3', 'Pool4'],
        'Sequence': [
            'CCTACACGACGCTCTTCCGATCTNNNNNNACATCAATGAAGACGGCATC',  # Forward primer
            'TCAGACGTGTGCTCTTCCGATCTNNNNNNATAGTTATGAAGGCTTCC',  # Reverse primer  
            'CCTACACGACGCTCTTCCGATCTNNNNNNACATCAATGAAGACGGCATC',  # Duplicate forward (same gene-specific part)
            'TCAGACGTGTGCTCTTCCGATCTNNNNNNATAGTTATCGGAATGCTCC',  # Different reverse
        ]
    }
    
    df = pd.DataFrame(test_data)
    
    with tempfile.NamedTemporaryFile(suffix='.xlsx', delete=False) as tmp_file:
        df.to_excel(tmp_file.name, index=False)
        excel_path = Path(tmp_file.name)
    
    try:
        yield excel_path
    finally:
        excel_path.unlink(missing_ok=True)


@pytest.fixture
def sample_primer_info() -> pd.DataFrame:
    """Create sample primer info DataFrame for testing."""
    return pd.DataFrame({
        'pool': ['Pool1', 'Pool2', 'Pool3', 'Pool4'],
        'primer_type': ['fwd', 'rev', 'fwd', 'rev'],
        'primer_seq_full': [
            'CCTACACGACGCTCTTCCGATCTNNNNNNACATCAATGAAGACGGCATC',
            'TCAGACGTGTGCTCTTCCGATCTNNNNNNATAGTTATGAAGGCTTCC',
            'CCTACACGACGCTCTTCCGATCTNNNNNNACATCATCCGATGCATTAG',
            'TCAGACGTGTGCTCTTCCGATCTNNNNNNATAGTTATCGGAATGCTCC',
        ],
        'primer_seq': [
            'ATGAAGACGGCATC',
            'ATGAAGGCTTCC',
            'ATCCGATGCATTAG',
            'ATCGGAATGCTCC',
        ]
    })


@pytest.fixture
def sample_locate_df() -> pd.DataFrame:
    """Create sample seqkit locate output DataFrame for testing."""
    return pd.DataFrame({
        'seqID': [
            'ENSMUST00000000001.4|ENSMUSG00000000001.4|ENSMUSG00000000001.4|3634||Gnai3|protein_coding',
            'ENSMUST00000000002.4|ENSMUSG00000000002.4|ENSMUSG00000000002.4|3635||Pbsn|protein_coding',
        ],
        'pattern': ['ATGAAGACGGCATC', 'ATGAAGGCTTCC'],
        'strand': ['+', '-'],
        'start': [100, 200],
        'end': [113, 211]
    })


@pytest.fixture
def sample_amplicon_df() -> pd.DataFrame:
    """Create sample amplicon DataFrame for testing."""
    return pd.DataFrame({
        'primer_seq_fwd': ['ATGAAGACGGCATC', 'ATCCGATGCATTAG'],
        'primer_seq_rev': ['ATGAAGGCTTCC', 'ATCGGAATGCTCC'],
        'transcript_id': ['ENSMUST00000000001.4', 'ENSMUST00000000002.4'],
        'gene_id': ['ENSMUSG00000000001.4', 'ENSMUSG00000000002.4'],
        'gene_symbol': ['Gnai3', 'Pbsn'],
        'start_up': [100, 150],
        'end_up': [113, 165],
        'start_down': [200, 250],
        'end_down': [211, 263],
        'pool_fwd': ['Pool1', 'Pool3'],
        'pool_rev': ['Pool2', 'Pool4'],
        'amplicon_length': [350, 360]
    })
