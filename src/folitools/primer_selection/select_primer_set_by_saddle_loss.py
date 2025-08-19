#!/usr/bin/env python
"""
Selects an optimal primer set from a list of candidates to minimize potential
primer-dimer formation.

This script implements a Simulated Annealing algorithm to solve a combinatorial
optimization problem. Given multiple candidate primer pairs for a list of genes,
it intelligently searches for the combination that results in the lowest
overall interaction score (SADDLE loss). The efficiency of the algorithm relies
on extensive pre-computation and a set of purpose-built data structures that allow
for rapid, incremental updates of the total loss score.

Algorithm Overview:
  1. Pre-computation: All expensive interaction calculations for every individual
     primer are performed once and stored in lookup tables (dictionaries).
  2. Initialization: An initial primer set is chosen (e.g., the first
     candidate for each gene) and its total SADDLE loss is calculated.
  3. Simulated Annealing: The script iteratively proposes a small change (swapping
     one primer for another) and probabilistically decides whether to accept the
     change based on the new score and the current "temperature". This allows
     the search to escape local minima.
  4. Refinement: A final, short "greedy" search is performed, only accepting
     changes that strictly improve the score to fine-tune the result.
  5. Output: The final, optimized primer set and the history of the SADDLE loss
     are saved to files.

Key Data Structures:
    The performance of this script hinges on the following data structures, which
    separate the static, pre-computed properties of a primer from the dynamic,
    system-wide state of the primer pool.

    primer_info_pool:
        - Type: dict[str, list[tuple[str, str, str]]]
        - Stores: The master database of all candidate primer pairs, loaded from
          the input file.
        - Structure: {gene_id: [(design_1, fwd_seq_1, rev_seq_1), ...]}
        - Update: Created once at the start by `load_primer_data` and is
          read-only thereafter.

    current_primer_used:
        - Type: list[int]
        - Stores: The core "state vector" of the simulation. It defines the
          exact combination of primers currently being evaluated.
        - Structure: An index-based list where `current_primer_used[i]` is the
          design index for the i-th gene.
        - Update: Initialized with a default set (e.g., all zeros). It is
          updated during the annealing loop whenever a new combination is
          accepted as the new "current best" state.

    dimer_length_gc:
        - Type: dict[str, float]
        - Stores: A fundamental scoring table that maps any short DNA sequence
          (a "tail") to a raw interaction score.
        - Structure: {tail_sequence: 2^(length + GC_count)}
        - Update: Created once during initialization. It is read-only.

    primer2selfloss:
        - Type: dict[str, float]
        - Stores: The pre-computed self-dimerization score for every unique
          primer sequence.
        - Structure: {primer_sequence: float_score}
        - Update: Created once during the pre-computation phase. It is
          read-only.

    primer2tail_score (The "Attacker"):
        - Type: dict[str, list[tuple[str, float]]]
        - Stores: A primer's "scoring checklist." It contains information about
          the primer's OWN 3' tails and their intrinsic binding scores. This
          is used to calculate how a primer interacts with the rest of the pool.
        - Structure: {primer_sequence: [(tail_1, score_1), (tail_2, score_2), ...]}
        - Update: Created once during pre-computation. It is read-only. Used by
          `calculate_saddle_by_hash`.

    primer2tailrc_weight (The "Defender"):
        - Type: dict[str, list[tuple[str, float]]]
        - Stores: A primer's "update signature." It describes the binding sites
          a primer presents to the pool, based on its REVERSE COMPLEMENT. This
          is used to modify the pool's overall state.
        - Structure: {primer_sequence: [(rc_tail_1, weight_1), ...]}
        - Update: Created once during pre-computation. It is read-only. Used by
          `update_tail_weight` to modify `current_tail2weight`.

    current_tail2weight:
        - Type: dict[str, float]
        - Stores: The "master compatibility ledger" for the entire primer set
          currently defined by `current_primer_used`. It aggregates the
          "defensive" interaction potential from all primers in the pool.
        - Structure: {tail_sequence: accumulated_weight}
        - Update: This is the primary dynamically updated data structure. It is
          initialized based on the starting set. During the annealing loop, it
          is updated incrementally by `update_tail_weight` each time a move
          is accepted, by subtracting the `primer2tailrc_weight` of the old
          primer and adding that of the new primer.
"""

from collections import defaultdict
import math
import random
from copy import deepcopy
from pathlib import Path
from functools import lru_cache

import pandas as pd

from .saddle_utils import choice_except


@lru_cache(maxsize=None)
def _tail_score(tail: str) -> float:
    """Calculate the score for a given tail sequence."""
    return 2 ** (len(tail) + count_gc(tail))


def init_hash_table(
    primer_pool: list[str], overlap_min: int, overlap_max: int, p3_dist_max: int
) -> dict[str, float]:
    """Allocate a tail-hash with zeros for both tails and RC tails.

    Args:
        primer_pool: Primer sequences (5'→3').
        overlap_min: Minimum overlap length to hash.
        overlap_max: Maximum overlap length to hash.
        p3_dist_max: Max offset from the 3' end to shift the tail.

    Returns:
        Mapping tail sequence -> 0.0 (ready to be incremented/decremented).
    """
    tail2weight_pool: dict[str, float] = {}
    for primer in primer_pool:
        rc = reverse_complement(primer)
        primer_len = len(primer)
        for overlap_len in range(overlap_min, overlap_max + 1):
            for p3_distance in range(p3_dist_max + 1):
                end = primer_len - p3_distance
                start = end - overlap_len
                if start < 0:
                    continue
                tail = primer[start:end]
                rc_tail = rc[p3_distance : p3_distance + overlap_len]
                tail2weight_pool[tail] = 0.0
                tail2weight_pool[rc_tail] = 0.0
    return tail2weight_pool


def reverse_complement(sequence: str) -> str:
    """Return the reverse-complement of a DNA sequence (ACGT only)."""
    pool = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(pool[b] for b in reversed(sequence))


def count_gc(seq: str) -> int:
    """Count the number of G/C bases in a sequence."""
    return sum(1 for base in seq if base in "GC")


def generate_primer_seq(
    genes: list[str],
    primer_info_pool: dict[str, list[tuple[str, str, str]]],
    design_idxs: list[int],
) -> list[str]:
    """Flatten the current selection into an alternating [F, R, F, R, ...] list.

    Args:
        genes: Gene IDs (order defines selection vector indexing).
        primer_info_pool: {gene -> [(design, fwd, rev), ...]}.
        design_idxs: For each gene index, which design index is chosen.

    Returns:
        List of sequences: [fwd_0, rev_0, fwd_1, rev_1, ...].
    """
    primer_seq_list: list[str] = []
    for gene, design_idx in zip(genes, design_idxs, strict=True):
        design, seq_f, seq_r = primer_info_pool[gene][design_idx]
        primer_seq_list.append(seq_f)
        primer_seq_list.append(seq_r)
    return primer_seq_list


def calc_tail_score(
    primer: str,
    overlap_min: int,
    overlap_max: int,
    p3_dist_max: int,
) -> list[tuple[str, float]]:
    """Build a compressed view of a primer's own 3' tails with scores.

    Args:
        primer_sequence: Primer (5'→3').
        dimer_length_gc: Tail -> 2**(len + GC).
        overlap_min: Min tail length.
        overlap_max: Max tail length.
        p3_dist_max: Max offset from 3' end.

    Returns:
        List of (tail, score) used for fast SADDLE accumulation.

    Examples:
        Using primer "ATGC" with overlap lengths 2-3 and p3 shifts up to 1.
        Assume ``dimer_length_gc[tail] = 2 ** (len(tail) + GC_count(tail))``.
        For the tails encountered here:
            - "GC"  -> 2**(2 + 2) = 16
            - "TG"  -> 2**(2 + 1) = 8
            - "TGC" -> 2**(3 + 2) = 32
            - "ATG" -> 2**(3 + 1) = 16

        >>> dimer_length_gc = {"GC": 16.0, "TG": 8.0, "TGC": 32.0, "ATG": 16.0}
        >>> pre_calculate_primer_self_hash(
        ...     "ATGC", dimer_length_gc, overlap_min=2, overlap_max=3, p3_dist_max=1
        ... )
        [('GC', 16.0), ('TG', 4.0), ('TGC', 32.0), ('ATG', 8.0)]
    """
    primer_len = len(primer)
    self_hash: list[tuple[str, float]] = []

    for overlap_len in range(overlap_min, overlap_max + 1):
        for p3_distance in range(p3_dist_max + 1):
            end = primer_len - p3_distance
            start = end - overlap_len
            if start < 0:
                continue
            tail = primer[start:end]
            # if tail in dimer_length_gc:
            score = _tail_score(tail) / (p3_distance + 1.0)
            self_hash.append((tail, score))
    return self_hash


def calc_tailrc_weight(
    primer: str, overlap_min: int, overlap_max: int, p3_dist_max: int
) -> list[tuple[str, float]]:
    """Compute RC-tail deltas this primer contributes to the global hash.

    Args:
        primer: Primer (5'→3').
        overlap_min: Min RC-tail length.
        overlap_max: Max RC-tail length.
        p3_dist_max: Max offset from 3' end.

    Returns:
        List of (rc_tail, weight) entries used to update the global hash.

    Examples:
        The reverse complement of "ATGC" is "GCAT". For overlap lengths 2-3
        and p3 shifts up to 1, the affected RC tails and their weights are:

        >>> calc_tailrc_weight(
        ...     "ATGC", overlap_min=2, overlap_max=3, p3_dist_max=1
        ... )
        [('GC', 1.0), ('CA', 0.5), ('GCA', 1.0), ('CAT', 0.5)]
    """
    primer_rc = reverse_complement(primer)
    affected_hash: list[tuple[str, float]] = []
    for overlap_len in range(overlap_min, overlap_max + 1):
        for p3_distance in range(p3_dist_max + 1):
            start = p3_distance
            end = start + overlap_len
            if end > len(primer_rc):
                continue
            rc_tail = primer_rc[start:end]
            weight = 1.0 / (p3_distance + 1)
            affected_hash.append((rc_tail, weight))
    return affected_hash


def calc_selfloss(
    primer: str, overlap_min: int, overlap_max: int, p3_dist_max: int
) -> float:
    """Compute self-dimer loss by matching 3' tails vs RC 5' windows.

    Args:
        primer: Primer (5'→3').
        overlap_min: Min overlap length.
        overlap_max: Max overlap length.
        p3_dist_max: Max 3' shift for both strands.

    Returns:
        Accumulated self-dimer loss.

    Examples:
        With primer "ATGC" (reverse complement "GCAT"), overlap lengths 2-3,
        and p3 shifts up to 1, the only exact match is tail "GC" at
        (p3_f=0, p3_r=0). Its contribution is
        ``2 ** (len('GC') + GC('GC')) / ((0 + 1) * (0 + 1)) = 2 ** (2 + 2) = 16.0``.

        >>> calc_selfloss(
        ...     "ATGC", overlap_min=2, overlap_max=3, p3_dist_max=1
        ... )
        16.0
    """
    primer_rc = reverse_complement(primer)
    primer_len = len(primer)
    total_loss = 0.0

    for overlap_len in range(overlap_min, overlap_max + 1):
        for p3_f in range(p3_dist_max + 1):
            end_f = primer_len - p3_f
            start_f = end_f - overlap_len
            if start_f < 0:
                continue
            tail_f = primer[start_f:end_f]

            for p3_r in range(p3_dist_max + 1):
                start_r = p3_r
                end_r = start_r + overlap_len
                if end_r > primer_len:
                    continue
                tail_r = primer_rc[start_r:end_r]

                if tail_f == tail_r:
                    gc = count_gc(tail_f)
                    loss = 2 ** (overlap_len + gc) / ((p3_f + 1) * (p3_r + 1))
                    total_loss += loss
    return total_loss


def update_tail_weight(
    tail2weight: dict[str, float] | None = None,
    *,
    old_primers: list[str],
    new_primers: list[str],
    primer2tailrc_weight: dict[str, list[tuple[str, float]]],
) -> dict:
    """Apply delta updates for removed then added primers.

    Args:
        tail2weight: Global hash (tail -> weight).
        old_primers: Primers to remove.
        new_primers: Primers to add.
        primer2tailrc_weight: Primer -> [(rc_tail, delta), ...].

    Returns:
        The same dict, modified in place and returned for convenience.
    """
    if tail2weight is None:
        tail2weight = defaultdict(float)
    for primer in old_primers:
        for tail, delta in primer2tailrc_weight[primer]:
            tail2weight[tail] -= delta
    for primer in new_primers:
        for tail, delta in primer2tailrc_weight[primer]:
            tail2weight[tail] += delta
    return tail2weight


def calc_saddle_by_hash(
    primers: list[str],
    tail2weight: dict[str, float],
    primer2tail_score: dict[str, list[tuple[str, float]]],
) -> float:
    """Accumulate SADDLE loss via precomputed tail scores and current hash.

    Args:
        primers: Flattened [F,R,...] list of chosen primers.
        tail2weight: Global hash mapping tail -> accumulated weight.
        primer2tail_score: Primer -> [(tail, score), ...].

    Returns:
        The interaction loss excluding the explicit self-loss term.
    """
    return sum(
        tail2weight[tail] * score
        for primer in primers
        for tail, score in primer2tail_score[primer]
    )


def calc_saddle_self_loss(
    old_primers: list[str],
    new_primers: list[str],
    primer2selfloss_table: dict[str, float],
) -> float:
    """Update the self-dimer loss component incrementally.

    Args:
        old_primers: Primers being removed.
        new_primers: Primers being added.
        primer2selfloss_table: Primer -> self-loss.

    Returns:
        New total self-loss.
    """
    new_loss = 0
    for p in old_primers:
        new_loss -= primer2selfloss_table[p]
    for p in new_primers:
        new_loss += primer2selfloss_table[p]
    return new_loss


def load_primer_data(
    input_: Path,
) -> tuple[dict[str, list[tuple[str, str, str]]], list[str]]:
    """Load candidate primers from TSV into grouped structures.

    Expected columns: gene, design, seq_f, seq_r (any header is coerced).

    Args:
        input_: TSV with candidate primer pairs.

    Returns:
        primer_info_pool: {gene -> [(design, seq_f, seq_r), ...]}
        primers:  Flattened list of all sequences for precomputation.
    """
    df = pd.read_csv(input_, sep="\t")
    df.columns = ["gene", "design", "seq_f", "seq_r"]
    print(f"Loaded {len(df)} primer records")

    primer_info_pool: dict[str, list[tuple[str, str, str]]] = {}
    primers: list[str] = []

    for gene, group in df.groupby("gene"):
        assert isinstance(gene, str)
        primer_info_pool[gene] = []
        for _, row in group.iterrows():
            design = row["design"]
            seq_f = row["seq_f"]
            seq_r = row["seq_r"]
            primer_info_pool[gene].append((design, seq_f, seq_r))
            primers.extend([seq_f, seq_r])

    print(
        f"Processed {len(primer_info_pool)} genes with {len(primers)} total sequences"
    )
    return primer_info_pool, primers


def saddle(
    *,
    input_: Path,
    output: Path,
    output_loss: Path,
    num_cycles_anneal: int = 50,
    random_seed: int = 42,
    overlap_min: int = 4,
    overlap_max: int = 12,
    p3_dist_max: int = 2,
    tolerance_factor: float = 1e6,
    reanneal_fraction: float = 0.5,
) -> int:
    """Run simulated annealing to choose a low-dimer primer set.

    Args:
        input_: TSV of candidate primers (gene, design, seq_f, seq_r).
        output: Output TSV with chosen primers.
        output_loss: Output TSV with loss trajectory (1 column, no header).
        num_cycles_anneal: Number of annealing iterations.
        random_seed: RNG seed.
        overlap_min: Minimum tail size for hashing (default: 4).
        overlap_max: Maximum tail size for hashing (default: 12).
        p3_dist_max: Max 3' offset considered (default: 2).
        tolerance_factor: SA acceptance scale (default: 1_000_000.0).
        reanneal_fraction: Fraction of cycles for greedy re-anneal (default: 0.5).

    Returns:
        Exit code 0 on success.
    """
    random.seed(random_seed)

    num_cycles_after = int(reanneal_fraction * num_cycles_anneal)
    print("Start!")

    primer_info_pool, primers = load_primer_data(input_)
    genes = list(primer_info_pool.keys())

    gene_total_design_list = [list(range(len(primer_info_pool[g]))) for g in genes]
    gene_idx_multi = [i for i, g in enumerate(genes) if len(primer_info_pool[g]) > 1]

    print("Initialize tail weights and scores and self losses...")

    primer2tail_score: dict[str, list[tuple[str, float]]] = {}
    primer2tailrc_weight: dict[str, list[tuple[str, float]]] = {}
    primer2selfloss: dict[str, float] = {}

    for p in primers:
        primer2tail_score[p] = calc_tail_score(p, overlap_min, overlap_max, p3_dist_max)
        primer2tailrc_weight[p] = calc_tailrc_weight(
            p, overlap_min, overlap_max, p3_dist_max
        )
        primer2selfloss[p] = calc_selfloss(p, overlap_min, overlap_max, p3_dist_max)

    # To start with, we select the first design for each gene.
    current_primer_used = [0] * len(genes)
    primers_init = generate_primer_seq(genes, primer_info_pool, current_primer_used)

    tail2weight_current = update_tail_weight(
        old_primers=[],
        new_primers=primers_init,
        primer2tailrc_weight=primer2tailrc_weight,
    )
    tail2weight_candidate = deepcopy(tail2weight_current)

    current_saddle_loss = calc_saddle_by_hash(
        primers=primers_init,
        tail2weight=tail2weight_current,
        primer2tail_score=primer2tail_score,
    ) + calc_saddle_self_loss(
        old_primers=[], new_primers=primers_init, primer2selfloss_table=primer2selfloss
    )
    list_saddle_loss = [current_saddle_loss]

    print("Start annealing primers...")
    for i in range(num_cycles_anneal):
        # Randomly choose a gene that has more than 1 design.
        gene_idx = random.choice(gene_idx_multi)
        gene_id = genes[gene_idx]

        cur_design = current_primer_used[gene_idx]
        cur_primer_seq = list(primer_info_pool[gene_id][cur_design][1:])

        new_design = choice_except(gene_total_design_list[gene_idx], cur_design)
        new_primer_seq = list(primer_info_pool[gene_id][new_design][1:])

        new_primer_used = current_primer_used[:]
        new_primer_used[gene_idx] = new_design

        update_tail_weight(
            tail2weight_candidate,
            old_primers=cur_primer_seq,
            new_primers=new_primer_seq,
            primer2tailrc_weight=primer2tailrc_weight,
        )

        candidate_loss = calc_saddle_by_hash(
            primers=generate_primer_seq(genes, primer_info_pool, new_primer_used),
            tail2weight=tail2weight_candidate,
            primer2tail_score=primer2tail_score,
        ) + calc_saddle_self_loss(
            old_primers=cur_primer_seq,
            new_primers=new_primer_seq,
            primer2selfloss_table=primer2selfloss,
        )

        loss_diff = candidate_loss - current_saddle_loss
        accept = False
        if loss_diff < 0:  # Always accept if loss is lower
            accept = True
        else:
            if random.random() < math.exp(-loss_diff / (tolerance_factor / (i + 1))):
                accept = True
        if accept:  # Accept, sync the change
            current_saddle_loss = candidate_loss
            current_primer_used = new_primer_used
            update_tail_weight(
                tail2weight_current,
                old_primers=cur_primer_seq,
                new_primers=new_primer_seq,
                primer2tailrc_weight=primer2tailrc_weight,
            )
            list_saddle_loss.append(current_saddle_loss)
        else:  # Do not accept, reverse the change
            update_tail_weight(
                tail2weight_candidate,
                old_primers=new_primer_seq,
                new_primers=cur_primer_seq,
                primer2tailrc_weight=primer2tailrc_weight,
            )
            list_saddle_loss.append(current_saddle_loss)

    print("Start re-annealing (greedy)...")
    for _ in range(int(num_cycles_after)):
        gene_idx = random.choice(gene_idx_multi)
        gene_id = genes[gene_idx]

        cur_design = current_primer_used[gene_idx]
        cur_primer_seq = list(primer_info_pool[gene_id][cur_design][1:])

        new_design = choice_except(gene_total_design_list[gene_idx], cur_design)
        new_primer_seq = list(primer_info_pool[gene_id][new_design][1:])

        new_primer_used = current_primer_used[:]
        new_primer_used[gene_idx] = new_design

        update_tail_weight(
            tail2weight_candidate,
            old_primers=cur_primer_seq,
            new_primers=new_primer_seq,
            primer2tailrc_weight=primer2tailrc_weight,
        )

        candidate_loss = calc_saddle_by_hash(
            primers=generate_primer_seq(genes, primer_info_pool, new_primer_used),
            tail2weight=tail2weight_candidate,
            primer2tail_score=primer2tail_score,
        ) + calc_saddle_self_loss(
            old_primers=cur_primer_seq,
            new_primers=new_primer_seq,
            primer2selfloss_table=primer2selfloss,
        )

        if candidate_loss < current_saddle_loss:
            current_saddle_loss = candidate_loss
            current_primer_used = new_primer_used
            update_tail_weight(
                tail2weight_current,
                old_primers=cur_primer_seq,
                new_primers=new_primer_seq,
                primer2tailrc_weight=primer2tailrc_weight,
            )
            list_saddle_loss.append(current_saddle_loss)
        else:
            update_tail_weight(
                tail2weight_candidate,
                old_primers=new_primer_seq,
                new_primers=cur_primer_seq,
                primer2tailrc_weight=primer2tailrc_weight,
            )
            list_saddle_loss.append(current_saddle_loss)
    print("Finished!")

    loss_df = pd.DataFrame({"saddle_loss": list_saddle_loss})
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
    print(f"  - SADDLE loss: {output_loss} ({len(loss_df)} values)")

    return 0
