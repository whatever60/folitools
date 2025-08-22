"""Tests for the eval_dimer module."""

import tempfile
from pathlib import Path
from typing import Generator

import numpy as np
import pandas as pd
import pytest

from folitools.primer_selection.eval_dimer import (
    _read_fasta,
    _pairwise_heterodimer_matrix,
    dimer_thermo_property,
)


@pytest.fixture
def sample_i5_fasta() -> Generator[Path, None, None]:
    """Create a temporary FASTA file with sample i5/forward primers."""
    fasta_content = """>i5_primer_001
ATGCGTACGGTA
>i5_primer_002
CGATCGATCGAT
>i5_primer_003
TTTTAAAACCCC
"""
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".fasta", delete=False
    ) as tmp_file:
        tmp_file.write(fasta_content)
        fasta_path = Path(tmp_file.name)

    try:
        yield fasta_path
    finally:
        fasta_path.unlink(missing_ok=True)


@pytest.fixture
def sample_i7_fasta() -> Generator[Path, None, None]:
    """Create a temporary FASTA file with sample i7/reverse primers."""
    fasta_content = """>i7_primer_001
GGGGTTTTAAAA
>i7_primer_002
ATCGATCGATCG
"""
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".fasta", delete=False
    ) as tmp_file:
        tmp_file.write(fasta_content)
        fasta_path = Path(tmp_file.name)

    try:
        yield fasta_path
    finally:
        fasta_path.unlink(missing_ok=True)


@pytest.fixture
def empty_fasta() -> Generator[Path, None, None]:
    """Create a temporary empty FASTA file."""
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".fasta", delete=False
    ) as tmp_file:
        tmp_file.write("")
        fasta_path = Path(tmp_file.name)

    try:
        yield fasta_path
    finally:
        fasta_path.unlink(missing_ok=True)


@pytest.fixture
def fasta_with_whitespace() -> Generator[Path, None, None]:
    """Create a FASTA file with sequences containing whitespace and newlines."""
    fasta_content = """>seq_with_spaces
ATG CGT ACG
 GTA
>seq_with_tabs
CGAT	CGAT
CGAT
"""
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".fasta", delete=False
    ) as tmp_file:
        tmp_file.write(fasta_content)
        fasta_path = Path(tmp_file.name)

    try:
        yield fasta_path
    finally:
        fasta_path.unlink(missing_ok=True)


class TestReadFasta:
    """Test the _read_fasta function."""

    def test_read_fasta_basic(self, sample_i5_fasta):
        """Test basic FASTA reading functionality."""
        names, seqs = _read_fasta(sample_i5_fasta)

        assert len(names) == 3
        assert len(seqs) == 3
        assert names == ["i5_primer_001", "i5_primer_002", "i5_primer_003"]
        assert seqs == ["ATGCGTACGGTA", "CGATCGATCGAT", "TTTTAAAACCCC"]

    def test_read_fasta_uppercase_conversion(self, fasta_with_whitespace):
        """Test that sequences are properly uppercased and whitespace is removed."""
        names, seqs = _read_fasta(fasta_with_whitespace)

        assert len(names) == 2
        assert len(seqs) == 2
        assert names == ["seq_with_spaces", "seq_with_tabs"]
        assert seqs == ["ATGCGTACGGTA", "CGATCGATCGAT"]
        # Check that whitespace has been removed
        assert " " not in seqs[0]
        assert "\t" not in seqs[1]
        assert "\n" not in seqs[0]
        assert "\n" not in seqs[1]

    def test_read_fasta_file_not_found(self):
        """Test FileNotFoundError when FASTA file doesn't exist."""
        nonexistent_path = Path("/nonexistent/path/file.fasta")

        with pytest.raises(FileNotFoundError, match="FASTA not found"):
            _read_fasta(nonexistent_path)

    def test_read_fasta_empty_file(self, empty_fasta):
        """Test ValueError when FASTA file is empty."""
        with pytest.raises(ValueError, match="No sequences found"):
            _read_fasta(empty_fasta)

    def test_read_fasta_pathlib_and_string_inputs(self, sample_i5_fasta):
        """Test that both Path objects and string paths work."""
        # Test with Path object
        names1, seqs1 = _read_fasta(sample_i5_fasta)

        # Test with string path
        names2, seqs2 = _read_fasta(str(sample_i5_fasta))

        assert names1 == names2
        assert seqs1 == seqs2


class TestPairwiseHeterodimerMatrix:
    """Test the _pairwise_heterodimer_matrix function."""

    def test_pairwise_heterodimer_matrix_shape(self):
        """Test that the output matrices have correct shape."""
        seqs_a = ["ATGCGTACGGTA", "CGATCGATCGAT"]
        seqs_b = ["GGGGTTTTAAAA", "ATCGATCGATCG", "TTTTAAAACCCC"]

        dg, tm = _pairwise_heterodimer_matrix(seqs_a, seqs_b)

        assert dg.shape == (2, 3)
        assert tm.shape == (2, 3)
        assert isinstance(dg, np.ndarray)
        assert isinstance(tm, np.ndarray)

    def test_pairwise_heterodimer_matrix_dtype(self):
        """Test that output matrices have correct data type."""
        seqs_a = ["ATGCGTACGGTA"]
        seqs_b = ["GGGGTTTTAAAA"]

        dg, tm = _pairwise_heterodimer_matrix(seqs_a, seqs_b)

        assert dg.dtype == float
        assert tm.dtype == float

    def test_pairwise_heterodimer_matrix_single_pair(self):
        """Test with a single sequence pair."""
        seqs_a = ["ATGCGTACGGTA"]
        seqs_b = ["GGGGTTTTAAAA"]

        dg, tm = _pairwise_heterodimer_matrix(seqs_a, seqs_b)

        assert dg.shape == (1, 1)
        assert tm.shape == (1, 1)
        # Values should be finite numbers
        assert np.isfinite(dg[0, 0])
        assert np.isfinite(tm[0, 0])

    def test_pairwise_heterodimer_matrix_empty_input(self):
        """Test with empty sequence lists."""
        seqs_a = []
        seqs_b = []

        dg, tm = _pairwise_heterodimer_matrix(seqs_a, seqs_b)

        assert dg.shape == (0, 0)
        assert tm.shape == (0, 0)


class TestDimerThermoProperty:
    """Test the dimer_thermo_property function."""

    def test_dimer_thermo_property_basic(self, sample_i5_fasta, sample_i7_fasta):
        """Test basic functionality of dimer thermodynamic calculation."""
        dg_df, tm_df = dimer_thermo_property(sample_i5_fasta, sample_i7_fasta)

        # Check DataFrame structure
        expected_size = 3 + 2  # 3 i5 primers + 2 i7 primers
        assert dg_df.shape == (expected_size, expected_size)
        assert tm_df.shape == (expected_size, expected_size)

        # Check index and column labels
        expected_labels = [
            "i5_primer_001_fwd",
            "i5_primer_002_fwd",
            "i5_primer_003_fwd",
            "i7_primer_001_rev",
            "i7_primer_002_rev",
        ]
        assert list(dg_df.index) == expected_labels
        assert list(dg_df.columns) == expected_labels
        assert list(tm_df.index) == expected_labels
        assert list(tm_df.columns) == expected_labels

        # Check that all values are finite numbers
        assert np.all(np.isfinite(dg_df.values))
        assert np.all(np.isfinite(tm_df.values))

    def test_dimer_thermo_property_custom_suffixes(
        self, sample_i5_fasta, sample_i7_fasta
    ):
        """Test with custom suffixes for primer names."""
        dg_df, tm_df = dimer_thermo_property(
            sample_i5_fasta, sample_i7_fasta, i5_suffix="forward", i7_suffix="reverse"
        )

        # Check that custom suffixes are used
        expected_labels = [
            "i5_primer_001_forward",
            "i5_primer_002_forward",
            "i5_primer_003_forward",
            "i7_primer_001_reverse",
            "i7_primer_002_reverse",
        ]
        assert list(dg_df.index) == expected_labels
        assert list(tm_df.index) == expected_labels

    def test_dimer_thermo_property_output_files(self, sample_i5_fasta, sample_i7_fasta):
        """Test CSV output functionality."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            dg_output = Path(tmp_dir) / "dg_output.csv"
            tm_output = Path(tmp_dir) / "tm_output.csv"

            dg_df, tm_df = dimer_thermo_property(
                sample_i5_fasta,
                sample_i7_fasta,
                output_dg_csv=dg_output,
                output_tm_csv=tm_output,
            )

            # Check that files were created
            assert dg_output.exists()
            assert tm_output.exists()

            # Check that files can be read back and match
            dg_loaded = pd.read_csv(dg_output, index_col=0)
            tm_loaded = pd.read_csv(tm_output, index_col=0)

            pd.testing.assert_frame_equal(dg_df, dg_loaded)
            pd.testing.assert_frame_equal(tm_df, tm_loaded)

    def test_dimer_thermo_property_output_directory_creation(
        self, sample_i5_fasta, sample_i7_fasta
    ):
        """Test that output directories are created if they don't exist."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            nested_dir = Path(tmp_dir) / "nested" / "directory"
            dg_output = nested_dir / "dg_output.csv"

            dg_df, tm_df = dimer_thermo_property(
                sample_i5_fasta, sample_i7_fasta, output_dg_csv=dg_output
            )

            # Check that nested directory was created
            assert nested_dir.exists()
            assert dg_output.exists()

    def test_dimer_thermo_property_file_not_found(self, sample_i5_fasta):
        """Test FileNotFoundError when input FASTA doesn't exist."""
        nonexistent_path = Path("/nonexistent/path/file.fasta")

        with pytest.raises(FileNotFoundError, match="FASTA not found"):
            dimer_thermo_property(sample_i5_fasta, nonexistent_path)

        with pytest.raises(FileNotFoundError, match="FASTA not found"):
            dimer_thermo_property(nonexistent_path, sample_i5_fasta)

    def test_dimer_thermo_property_empty_fasta(self, sample_i5_fasta, empty_fasta):
        """Test ValueError when one of the FASTA files is empty."""
        with pytest.raises(ValueError, match="No sequences found"):
            dimer_thermo_property(sample_i5_fasta, empty_fasta)

        with pytest.raises(ValueError, match="No sequences found"):
            dimer_thermo_property(empty_fasta, sample_i5_fasta)

    def test_dimer_thermo_property_matrix_blocks(
        self, sample_i5_fasta, sample_i7_fasta
    ):
        """Test that the combined matrix has correct block structure."""
        dg_df, tm_df = dimer_thermo_property(sample_i5_fasta, sample_i7_fasta)

        n_i5 = 3
        n_i7 = 2

        # Top-left block: i5 x i5
        i5_block_dg = dg_df.iloc[:n_i5, :n_i5]
        i5_block_tm = tm_df.iloc[:n_i5, :n_i5]
        assert i5_block_dg.shape == (n_i5, n_i5)
        assert i5_block_tm.shape == (n_i5, n_i5)

        # Bottom-right block: i7 x i7
        i7_block_dg = dg_df.iloc[n_i5:, n_i5:]
        i7_block_tm = tm_df.iloc[n_i5:, n_i5:]
        assert i7_block_dg.shape == (n_i7, n_i7)
        assert i7_block_tm.shape == (n_i7, n_i7)

        # Top-right block: i5 x i7
        i5i7_block_dg = dg_df.iloc[:n_i5, n_i5:]
        i5i7_block_tm = tm_df.iloc[:n_i5, n_i5:]
        assert i5i7_block_dg.shape == (n_i5, n_i7)
        assert i5i7_block_tm.shape == (n_i5, n_i7)

        # Bottom-left block: i7 x i5
        i7i5_block_dg = dg_df.iloc[n_i5:, :n_i5]
        i7i5_block_tm = tm_df.iloc[n_i5:, :n_i5]
        assert i7i5_block_dg.shape == (n_i7, n_i5)
        assert i7i5_block_tm.shape == (n_i7, n_i5)

    def test_dimer_thermo_property_string_paths(self, sample_i5_fasta, sample_i7_fasta):
        """Test that string paths work as well as Path objects."""
        # Test with Path objects
        dg_df1, tm_df1 = dimer_thermo_property(sample_i5_fasta, sample_i7_fasta)

        # Test with string paths
        dg_df2, tm_df2 = dimer_thermo_property(
            str(sample_i5_fasta), str(sample_i7_fasta)
        )

        pd.testing.assert_frame_equal(dg_df1, dg_df2)
        pd.testing.assert_frame_equal(tm_df1, tm_df2)

    def test_dimer_thermo_property_single_primer_files(self):
        """Test with FASTA files containing only one primer each."""
        single_i5_content = """>single_i5
ATGCGTACGGTA
"""
        single_i7_content = """>single_i7
GGGGTTTTAAAA
"""

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        ) as tmp_i5:
            tmp_i5.write(single_i5_content)
            i5_path = Path(tmp_i5.name)

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        ) as tmp_i7:
            tmp_i7.write(single_i7_content)
            i7_path = Path(tmp_i7.name)

        try:
            dg_df, tm_df = dimer_thermo_property(i5_path, i7_path)

            # Should have 2x2 matrix (1 i5 + 1 i7)
            assert dg_df.shape == (2, 2)
            assert tm_df.shape == (2, 2)

            expected_labels = ["single_i5_fwd", "single_i7_rev"]
            assert list(dg_df.index) == expected_labels
            assert list(dg_df.columns) == expected_labels

        finally:
            i5_path.unlink(missing_ok=True)
            i7_path.unlink(missing_ok=True)
