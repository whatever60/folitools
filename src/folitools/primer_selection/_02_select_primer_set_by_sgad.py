#!/usr/bin/env python3
"""Simulated annealing primer selection using SGAD dimer scoring.

Scoring convention is symmetric and unordered at the primer-instance level:
- Each unordered primer-instance pair is counted once.
- Self pairs are included once.
- (primer_a, primer_b) and (primer_b, primer_a) are duplicates.
"""

from __future__ import annotations

import math
import random
from pathlib import Path

import pandas as pd
from sgad.rust.pairwise import make_rust_score_scaler
from sgad.rust.pairwise import needleman_wunsch_batch as rust_needleman_wunsch_batch
from tqdm.auto import trange

from ._02_select_primer_set_by_saddle_loss import (
    generate_primer_seq,
    load_primer_data,
    read_fasta,
)
from .saddle_utils import choice_except


_DNA_COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def _score_matrix(
    *, at_match: float, gc_match: float, mismatch: float
) -> dict[str, dict[str, float]]:
    """Build a base-aware substitution matrix for SGAD alignment."""
    return {
        "A": {"A": at_match, "C": mismatch, "G": mismatch, "T": mismatch},
        "C": {"A": mismatch, "C": gc_match, "G": mismatch, "T": mismatch},
        "G": {"A": mismatch, "C": mismatch, "G": gc_match, "T": mismatch},
        "T": {"A": mismatch, "C": mismatch, "G": mismatch, "T": at_match},
    }


def _replace_once(
    pool: list[str], old_primers: list[str], new_primers: list[str]
) -> list[str]:
    """Return a copy of pool where each old primer instance is removed once."""
    out = pool.copy()
    for primer in old_primers:
        out.remove(primer)
    out.extend(new_primers)
    return out


def _score_seq_pairs(
    seq_pairs: list[tuple[str, str]],
    *,
    matrix: dict[str, dict[str, float]],
    gap_open: float,
    gap_extend: float,
    score_scale_fn,
) -> float:
    """Score all sequence pairs in one Rust batch call."""
    if not seq_pairs:
        return 0.0

    scored = rust_needleman_wunsch_batch(
        [(seq1, seq2.translate(_DNA_COMP)[::-1]) for seq1, seq2 in seq_pairs],
        score_matrix=matrix,
        gap_open=gap_open,
        gap_extend=gap_extend,
        seq1_left_free=False,
        seq1_right_free=True,
        seq2_left_free=True,
        seq2_right_free=False,
        score_scale_fn=score_scale_fn,
    )
    return float(sum(float(score) for _, _, score in scored))


def _score_pool_full(
    pool: list[str],
    *,
    matrix: dict[str, dict[str, float]],
    gap_open: float,
    gap_extend: float,
    score_scale_fn,
) -> float:
    """Compute full SGAD pool score from all unordered primer-instance pairs."""
    pairs: list[tuple[str, str]] = []
    for i, seq1 in enumerate(pool):
        for seq2 in pool[i:]:
            pairs.append((seq1, seq2))
    return _score_seq_pairs(
        pairs,
        matrix=matrix,
        gap_open=gap_open,
        gap_extend=gap_extend,
        score_scale_fn=score_scale_fn,
    )


def _indices_for_remove_once(pool: list[str], to_remove: list[str]) -> list[int]:
    """Return indices removed by list.remove-style one-by-one primer deletion."""
    removed_indices: list[int] = []
    used_indices: set[int] = set()
    for seq in to_remove:
        for idx, pool_seq in enumerate(pool):
            if pool_seq == seq and idx not in used_indices:
                removed_indices.append(idx)
                used_indices.add(idx)
                break
    return removed_indices


def _score_affected_indices(
    pool: list[str],
    affected_indices: list[int],
    *,
    matrix: dict[str, dict[str, float]],
    gap_open: float,
    gap_extend: float,
    score_scale_fn,
) -> float:
    """Score all unordered pairs where at least one endpoint index is affected."""
    seen_pairs: set[tuple[int, int]] = set()
    seq_pairs: list[tuple[str, str]] = []

    for affected_idx in affected_indices:
        for idx, seq in enumerate(pool):
            pair_idx = (idx, affected_idx) if idx <= affected_idx else (affected_idx, idx)
            if pair_idx in seen_pairs:
                continue
            seen_pairs.add(pair_idx)
            seq_pairs.append((pool[pair_idx[0]], pool[pair_idx[1]]))

    return _score_seq_pairs(
        seq_pairs,
        matrix=matrix,
        gap_open=gap_open,
        gap_extend=gap_extend,
        score_scale_fn=score_scale_fn,
    )


def sgad(
    *,
    input_: Path,
    output: Path,
    output_loss: Path,
    num_cycles_anneal: int = 50,
    random_seed: int = 42,
    at_match: float,
    gc_match: float,
    mismatch: float,
    gap_open: float,
    gap_extend: float,
    decay_exponent: float,
    temperature: float,
    tolerance_factor: float = 1e6,
    reanneal_fraction: float = 0.5,
    skip_initial_pairwise_loss: bool = True,
    background_fasta: Path | None = None,
) -> pd.DataFrame:
    """Run SGAD-based simulated annealing to select one design per gene.

    Lower objective is better.

    If skip_initial_pairwise_loss is True, the initial score is set to zero and
    optimization proceeds by delta updates only.
    """
    random.seed(random_seed)
    num_cycles_after = int(reanneal_fraction * num_cycles_anneal)

    print("Start SGAD optimization!")
    primer_info_pool, _ = load_primer_data(input_)

    background_primers: list[str] = []
    if background_fasta is not None:
        background_primers = list(read_fasta(background_fasta).values())

    genes = list(primer_info_pool.keys())
    gene_total_design_list = [list(range(len(primer_info_pool[g]))) for g in genes]
    gene_idx_multi = [i for i, g in enumerate(genes) if len(primer_info_pool[g]) > 1]

    current_primer_used = [0] * len(genes)
    selected_primers = generate_primer_seq(genes, primer_info_pool, current_primer_used)
    current_pool = selected_primers + background_primers

    matrix = _score_matrix(at_match=at_match, gc_match=gc_match, mismatch=mismatch)
    score_scale_fn = make_rust_score_scaler(
        decay_exponent=decay_exponent,
        temperature=temperature,
    )

    if skip_initial_pairwise_loss:
        current_sgad_loss = 0.0
    else:
        current_sgad_loss = _score_pool_full(
            current_pool,
            matrix=matrix,
            gap_open=gap_open,
            gap_extend=gap_extend,
            score_scale_fn=score_scale_fn,
        )

    list_sgad_loss = [current_sgad_loss]

    print("Start simulated annealing (SGAD objective)...")
    prog_bar = trange(num_cycles_anneal + num_cycles_after)
    for i in prog_bar:
        gene_idx = random.choice(gene_idx_multi)
        gene_id = genes[gene_idx]

        cur_design = current_primer_used[gene_idx]
        cur_primer_seq = list(primer_info_pool[gene_id][cur_design][1:])

        new_design = choice_except(gene_total_design_list[gene_idx], cur_design)
        new_primer_seq = list(primer_info_pool[gene_id][new_design][1:])

        candidate_primer_used = current_primer_used[:]
        candidate_primer_used[gene_idx] = new_design
        candidate_pool = _replace_once(current_pool, cur_primer_seq, new_primer_seq)

        removed_indices = _indices_for_remove_once(current_pool, cur_primer_seq)
        added_indices = list(range(len(candidate_pool) - len(new_primer_seq), len(candidate_pool)))

        drop_score = _score_affected_indices(
            current_pool,
            removed_indices,
            matrix=matrix,
            gap_open=gap_open,
            gap_extend=gap_extend,
            score_scale_fn=score_scale_fn,
        )
        add_score = _score_affected_indices(
            candidate_pool,
            added_indices,
            matrix=matrix,
            gap_open=gap_open,
            gap_extend=gap_extend,
            score_scale_fn=score_scale_fn,
        )

        candidate_loss = current_sgad_loss - drop_score + add_score

        loss_diff = candidate_loss - current_sgad_loss
        if loss_diff < 0 or (
            i < num_cycles_anneal
            and random.random() < math.exp(-loss_diff / (tolerance_factor / (i + 1)))
        ):
            current_sgad_loss = candidate_loss
            current_primer_used = candidate_primer_used
            current_pool = candidate_pool

        prog_bar.set_postfix_str(f"SGAD loss: {current_sgad_loss:.2f}")
        list_sgad_loss.append(current_sgad_loss)

    print("Finished SGAD optimization!")

    loss_df = pd.DataFrame({"sgad_loss": list_sgad_loss})
    loss_df.to_csv(output_loss, sep="\t", index=False, header=False)

    rows = []
    for i, gene in enumerate(genes):
        design_idx = current_primer_used[i]
        design, seq_f, seq_r = primer_info_pool[gene][design_idx]
        rows.append(
            {
                "transcriptID": gene,
                "primerIndex": design,
                "sequenceLeft": seq_f,
                "sequenceRight": seq_r,
            }
        )

    results_df = pd.DataFrame(rows)
    results_df.to_csv(output, sep="\t", index=False)

    print("Results saved:")
    print(f"  - Final primers: {output} ({len(results_df)} records)")
    print(f"  - SGAD loss: {output_loss} ({len(loss_df)} values)")

    return results_df
