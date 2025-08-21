from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
import primer3
from tqdm.auto import tqdm


def _read_fasta(path: str | Path) -> tuple[list[str], list[str]]:
    """Read sequences and names from a FASTA file.

    Args:
        path: Path to a FASTA file.

    Returns:
        A pair (names, seqs) where each sequence is uppercased and stripped of whitespace.

    Raises:
        FileNotFoundError: If the path does not exist.
        ValueError: If the file contains no sequences.
    """
    fasta_path = Path(path)
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA not found: {fasta_path}")

    names: list[str] = []
    seqs: list[str] = []
    for record in SeqIO.parse(str(fasta_path), "fasta"):
        names.append(str(record.id))
        seqs.append(
            str(record.seq).upper().replace(" ", "").replace("\n", "").replace("\r", "")
        )

    if not names:
        raise ValueError(f"No sequences found in {fasta_path}")

    return names, seqs


def _pairwise_heterodimer_matrix(
    seqs_a: list[str], seqs_b: list[str]
) -> tuple[np.ndarray, np.ndarray]:
    """Compute ΔG and Tm matrices for all pairs (a in seqs_a, b in seqs_b).

    Args:
        seqs_a: First sequence list (treated as 5'→3').
        seqs_b: Second sequence list (treated as 5'→3').

    Returns:
        (dg, tm) where each is a 2D NumPy array of shape (len(seqs_a), len(seqs_b)).
    """
    n_a = len(seqs_a)
    n_b = len(seqs_b)
    dg = np.empty((n_a, n_b), dtype=float)
    tm = np.empty((n_a, n_b), dtype=float)

    # Nested loops are deliberate: primer3 works per-pair; vectorization isn't applicable.
    for i, sa in enumerate(tqdm(seqs_a)):
        for j, sb in enumerate(seqs_b):
            res = primer3.bindings.calc_heterodimer(sa, sb)
            dg[i, j] = res.dg
            tm[i, j] = res.tm

    return dg, tm


def dimer_thermo_property(
    path_i5_fasta: str | Path,
    path_i7_fasta: str | Path,
    *,
    i5_suffix: str = "fwd",
    i7_suffix: str = "rev",
    output_dg_csv: str | Path | None = None,
    output_tm_csv: str | Path | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compute primer-dimer thermodynamics for two primer sets (i5 and i7) and return combined matrices.

    This function assumes:
    - The first FASTA (``path_i5_fasta``) contains your i5 / forward primers.
    - The second FASTA (``path_i7_fasta``) contains your i7 / reverse primers.

    Orientation matters for primer-dimer calculations in `primer3`, so the resulting
    submatrices are generally **not symmetric**. We therefore compute:
      • i5xi5  (top-left block)
      • i7xi7  (bottom-right block)
      • i5xi7  (top-right block)
      • i7xi5  (bottom-left block)

    The combined matrices have index/columns:
        [<i5_name>_{i5_suffix}, ..., <i7_name>_{i7_suffix}, ...]

    Args:
        path_i5_fasta: FASTA path for i5 (forward) primers.
        path_i7_fasta: FASTA path for i7 (reverse) primers.
        i5_suffix: Suffix appended to i5 primer names in the combined matrices (default: "fwd").
        i7_suffix: Suffix appended to i7 primer names in the combined matrices (default: "rev").
        output_dg_csv: Optional CSV path to save the combined ΔG matrix.
        output_tm_csv: Optional CSV path to save the combined Tm matrix.

    Returns:
        A tuple (dg_combined_df, tm_combined_df), each a square pandas DataFrame whose rows/columns
        are the concatenated primer names from i5 then i7 with the chosen suffixes.

    Raises:
        FileNotFoundError: If any input FASTA path does not exist.
        ValueError: If a FASTA file is empty (no sequences).
    """
    # Read inputs
    i5_names, i5_seqs = _read_fasta(path_i5_fasta)
    i7_names, i7_seqs = _read_fasta(path_i7_fasta)

    n_i5 = len(i5_seqs)
    n_i7 = len(i7_seqs)

    # Compute all four blocks
    dg_i5, tm_i5 = _pairwise_heterodimer_matrix(i5_seqs, i5_seqs)  # top-left
    dg_i7, tm_i7 = _pairwise_heterodimer_matrix(i7_seqs, i7_seqs)  # bottom-right
    dg_i5i7, tm_i5i7 = _pairwise_heterodimer_matrix(i5_seqs, i7_seqs)  # top-right
    dg_i7i5, tm_i7i5 = _pairwise_heterodimer_matrix(i7_seqs, i5_seqs)  # bottom-left

    # Assemble combined matrices
    combined_size = n_i5 + n_i7
    dg_combined = np.zeros((combined_size, combined_size), dtype=float)
    tm_combined = np.zeros((combined_size, combined_size), dtype=float)

    # Block placement
    dg_combined[:n_i5, :n_i5] = dg_i5
    dg_combined[n_i5:, n_i5:] = dg_i7
    dg_combined[:n_i5, n_i5:] = dg_i5i7
    dg_combined[n_i5:, :n_i5] = dg_i7i5

    tm_combined[:n_i5, :n_i5] = tm_i5
    tm_combined[n_i5:, n_i5:] = tm_i7
    tm_combined[:n_i5, n_i5:] = tm_i5i7
    tm_combined[n_i5:, :n_i5] = tm_i7i5

    # Build labels: i5 first, then i7 (to match the block layout)
    idx_labels = [f"{name}_{i5_suffix}" for name in i5_names] + [
        f"{name}_{i7_suffix}" for name in i7_names
    ]

    dg_combined_df = pd.DataFrame(dg_combined, index=idx_labels, columns=idx_labels)
    tm_combined_df = pd.DataFrame(tm_combined, index=idx_labels, columns=idx_labels)

    # Optional saves
    if output_dg_csv is not None:
        Path(output_dg_csv).parent.mkdir(parents=True, exist_ok=True)
        dg_combined_df.to_csv(output_dg_csv)
    if output_tm_csv is not None:
        Path(output_tm_csv).parent.mkdir(parents=True, exist_ok=True)
        tm_combined_df.to_csv(output_tm_csv)

    return dg_combined_df, tm_combined_df
