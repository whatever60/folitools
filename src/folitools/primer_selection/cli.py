#!/usr/bin/env python3
"""
folitools.primer_selection CLI (cyclopts)

Subcommands
-----------
subset
    Subset packaged Parquet inputs by species + amplicon size range (TSV outputs).
saddle
    Select an optimal primer set using SADDLE (simulated annealing).
product
    Extract amplicon regions into FASTA using packaged references (by --species).

workflow
    Run subset → saddle → product end-to-end.

Examples
--------
folitools-primer subset \
  --species mouse \
  --amplicon-size-range 320 380 \
  --gene-table-file ./genes.tsv \
  --output-dir ./out

folitools-primer saddle \
  --input-file ./out/candidate_primer.320_380.primer_sequence.tsv \
  --output-final ./out/output_selected.tsv \
  --output-loss ./out/output_select.loss.txt \
  --num-cycles-anneal 50 \
  --random-seed 42 \
  --background-fasta ./background_primers.fasta

folitools-primer product \
  --input-selected ./out/output_selected.tsv \
  --primer-info ./out/candidate_primer.320_380.primer_info.tsv \
  --species mouse \
  --output-fasta ./out/output_selected.region.fasta

folitools-primer workflow \
  --input ./genes.tsv \
  --species mouse \
  --amplicon-size-range 320 380 \
  --num-cycles-anneal 50 \
  --random-seed 42 \
  --output-dir ./out \
  --background-fasta ./background_primers.fasta
"""

from pathlib import Path
from typing import Literal

from cyclopts import App

from ._01_primer_selection import subset as _subset
from ._02_select_primer_set_by_saddle_loss import saddle as _saddle
from ._03_extract_region_sequence import product as _product
from ._04_make_excel import summary as _summary


app = App(name="folitools-primer")


@app.command
def subset(
    *,
    input_: Path,
    species: str,
    amplicon_size_range: tuple[int, int],
    output_primer_sequence: Path,
    output_primer_info: Path,
) -> int:
    """Run the subset stage (packaged Parquet → filtered TSVs).

    Args:
        input_: TSV with columns: gene, group; optional primer_fwd, primer_rev.
        species: Species key (e.g., "mouse") mapping to packaged data.
        amplicon_size_range: Two integers MIN MAX. Example: `--amplicon-size-range 320 380`.
        output_primer_sequence: Output path for primer sequence TSV (must end with .tsv/.txt/.tsv.gz/.txt.gz).
        output_primer_info: Output path for primer info TSV (must end with .tsv/.txt/.tsv.gz/.txt.gz).

    Returns:
        Process exit code.
    """
    return _subset(
        species=species,
        amplicon_size_range=amplicon_size_range,
        gene_table_file=input_,
        output_primer_sequence=output_primer_sequence,
        output_primer_info=output_primer_info,
    )


@app.command
def saddle(
    *,
    input_: Path,
    output: Path,
    output_loss: Path,
    num_cycles_anneal: int = 50,
    random_seed: int = 42,
    background_fasta: Path | None = None,
) -> int:
    """Run SADDLE to choose a low-interaction primer set.

    Args:
        input_: TSV of candidate primers (gene, design, seq_f, seq_r).
        output: TSV path to write final chosen primers.
        output_loss: TSV path to write loss trajectory (no header).
        num_cycles_anneal: Annealing iterations.
        random_seed: RNG seed for reproducibility.
        background_fasta: Optional FASTA file with additional primers to include in the pool.

    Returns:
        Process exit code.
    """
    return _saddle(
        input_=input_,
        output=output,
        output_loss=output_loss,
        num_cycles_anneal=num_cycles_anneal,
        random_seed=random_seed,
        background_fasta=background_fasta,
    )


@app.command
def product(
    *,
    input_: Path,
    primer_info: Path,
    output_fasta: Path,
    species: Literal["mouse", "human"] | None = None,
    reference: Path | None = None,
) -> int:
    """Extract amplicon regions to FASTA from selected primer set.

    Args:
        input_: TSV from SADDLE stage (selected pairs).
        primer_info: TSV from subset stage (primer metadata).
        output_fasta: Output FASTA path for extracted regions.
        species: If provided (and `reference` not set), use packaged reference for species.
        reference: Optional explicit FASTA path (overrides `species`).

    Returns:
        Process exit code.
    """

    return _product(
        selected_tsv=input_,
        primer_info_tsv=primer_info,
        output_fasta=output_fasta,
        species=species,
        reference=reference,
    )


@app.command
def summary(
    *, input_: Path, primer_selection: Path, primer_info: Path, output: Path
) -> int:
    """Generate an Excel file for primer ordering.

    Args:
        input_: Path to gene info TSV file.
        primer_selection: Path to selected primers TSV file.
        primer_info: Path to candidate primer info TSV file.
        output: Path to output Excel file (default: primer_to_order.xlsx).

    Returns:
        Process exit code.
    """
    _summary(
        str(input_),
        str(primer_selection),
        str(primer_info),
        str(output),
    )
    return 0


@app.command
def workflow(
    *,
    input_: Path,  # gene table TSV
    species: Literal["mouse", "human"],
    reference: Path | None = None,
    amplicon_size_range: tuple[int, int],
    output_dir: Path,
    num_cycles_anneal: int = 50,
    random_seed: int = 42,
    background_fasta: Path | None = None,
) -> int:
    """Run the full pipeline: subset → saddle → product.

    Args:
        input_: Gene table TSV (columns: gene, group; optional primer_fwd, primer_rev).
        species: Species key for packaged data and reference (e.g., "mouse").
        amplicon_size_range: Two integers MIN MAX (e.g., 320 380).
        output_dir: Directory to write all outputs.
        num_cycles_anneal: SADDLE iterations (default: 50).
        random_seed: RNG seed (default: 42).
        background_fasta: Optional FASTA file with additional primers to include in the pool.

    Returns:
        Process exit code.
    """
    amin, amax = amplicon_size_range
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1) subset
    primer_seq = output_dir / "primer_sequence.tsv"
    primer_info = output_dir / "primer_info.tsv"
    rc1 = _subset(
        species=species,
        amplicon_size_range=amplicon_size_range,
        gene_table_file=input_,
        output_primer_sequence=primer_seq,
        output_primer_info=primer_info,
    )
    if rc1 != 0:
        return rc1

    # 2) saddle
    selected_tsv = output_dir / "primer_sequence_selected.tsv"
    loss_tsv = output_dir / "primer_sequence_selected_loss.txt"
    rc2 = saddle(
        input_=primer_seq,
        output=selected_tsv,
        output_loss=loss_tsv,
        num_cycles_anneal=num_cycles_anneal,
        random_seed=random_seed,
        background_fasta=background_fasta,
    )
    if rc2 != 0:
        return rc2

    # 3) product
    region_fa = output_dir / "amplicons.fasta"
    rc3 = _product(
        selected_tsv=selected_tsv,
        primer_info_tsv=primer_info,
        output_fasta=region_fa,
        species=species,
        reference=reference,
    )
    if rc3 != 0:
        return rc3

    # 4) summary (Excel output)
    excel_out = output_dir / "primer_to_order.xlsx"
    _summary(
        str(input_),
        str(selected_tsv),
        str(primer_info),
        str(excel_out),
    )
    return 0


def main() -> int:
    """CLI entrypoint."""
    try:
        app()
        return 0
    except SystemExit as e:
        return int(e.code or 0)


if __name__ == "__main__":
    raise SystemExit(main())
