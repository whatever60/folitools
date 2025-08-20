"""Shared utilities for primer selection modules."""

from __future__ import annotations

from pathlib import Path
from typing import Literal

from .data_dir_resolve import _data_dir_for_species

__all__ = ["get_prefixes", "resolve_reference_path", "SPECIES_FASTA_NAME"]

# Map species to packaged transcript FASTA filenames
SPECIES_FASTA_NAME: dict[str, str] = {
    "mouse": "Mus_musculus.GRCm39.transcript.fasta.gz",
    "human": "Homo_sapiens.GRCh38.transcript.fasta.gz",
}


def get_prefixes(has_linker: bool) -> tuple[str, str]:
    """Return (forward_prefix, reverse_prefix) depending on linker usage.

    If `has_linker` is True, include the 6bp linker (ACATCA/ATAGTT) after NNNNNN.

    Args:
        has_linker: Whether to include the 6bp linker sequences.

    Returns:
        Tuple of (forward_prefix, reverse_prefix) strings.
    """
    if has_linker:
        fwd = "CCTACACGACGCTCTTCCGATCTNNNNNNACATCA"
        rev = "TCAGACGTGTGCTCTTCCGATCTNNNNNNATAGTT"
    else:
        fwd = "CCTACACGACGCTCTTCCGATCTNNNNNN"
        rev = "TCAGACGTGTGCTCTTCCGATCTNNNNNN"
    return fwd, rev


def resolve_reference_path(
    txome_fasta: Path | None, species: Literal["mouse", "human"] | None
) -> Path:
    """Resolve transcriptome FASTA path from species or direct path.

    Args:
        txome_fasta: Path to the transcript FASTA (gz ok). If provided, overrides `species`.
        species: One of {"mouse", "human"} to resolve a packaged FASTA via `_data_dir_for_species`.

    Returns:
        Absolute :class:`Path` to the transcript FASTA.

    Raises:
        ValueError: If neither `species` nor `txome_fasta` is provided.
        FileNotFoundError: If the resolved FASTA path does not exist.
    """
    if species is not None:
        species_dir = _data_dir_for_species(species)
        ref_path = species_dir / SPECIES_FASTA_NAME[species]
    elif txome_fasta is not None:
        ref_path = Path(txome_fasta)
    else:
        raise ValueError("Either species or txome_fasta must be provided.")

    if not ref_path.exists():
        raise FileNotFoundError(f"Transcript FASTA not found: {ref_path}")
    return ref_path.resolve()
