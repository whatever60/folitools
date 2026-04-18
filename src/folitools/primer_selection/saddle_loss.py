"""
saddle_loss.py

Standalone SADDLE loss calculation for primer pools.

This module provides a simple interface to calculate the total SADDLE loss
for a given list of primer sequences without requiring pre-computation or
optimization infrastructure.
"""

from collections import defaultdict

from ._02_select_primer_set_by_saddle_loss import (
    calc_tail_score,
    calc_tailrc_weight,
    calc_selfloss,
)


def saddle_loss_naive(
    primers: list[str],
    *,
    overlap_min: int = 4,
    overlap_max: int = 12,
    p3_dist_max: int = 2,
) -> float:
    """Calculate SADDLE loss for a primer pool from scratch.

    This is a "naive" implementation that computes the full SADDLE loss without
    pre-computation or caching. It's useful for:
    - Understanding how the loss is calculated
    - Testing and validation
    - One-off calculations on small primer sets

    For optimization over many iterations, use the incremental update approach
    in `_02_select_primer_set_by_saddle_loss.py` instead.

    The SADDLE loss has two components:
    1. Interaction loss: How strongly primers bind to each other based on
       3' tail complementarity.
    2. Self loss: Self-dimerization potential of individual primers.

    Args:
        primers: List of primer sequences (5'->3', A/C/G/T only).
        overlap_min: Minimum tail length for scoring (default: 4).
        overlap_max: Maximum tail length for scoring (default: 12).
        p3_dist_max: Maximum offset from 3' end (default: 2).

    Returns:
        Total SADDLE loss (float).

    Examples:
        >>> primers = ["ATGCGATCGATCG", "CGATCGATCGCAT"]
        >>> loss = saddle_loss_naive(primers)
        >>> print(f"Total loss: {loss:.2f}")

        # With custom parameters
        >>> loss = saddle_loss_naive(
        ...     primers,
        ...     overlap_min=5,
        ...     overlap_max=10,
        ...     p3_dist_max=1
        ... )
    """
    if not primers:
        return 0.0

    # Step 1: Pre-compute properties for each unique primer
    primer2tail_score: dict[str, list[tuple[str, float]]] = {}
    primer2tailrc_weight: dict[str, list[tuple[str, float]]] = {}
    primer2selfloss: dict[str, float] = {}

    for p in primers:
        if p not in primer2tail_score:
            primer2tail_score[p] = calc_tail_score(p, overlap_min, overlap_max, p3_dist_max)
            primer2tailrc_weight[p] = calc_tailrc_weight(p, overlap_min, overlap_max, p3_dist_max)
            primer2selfloss[p] = calc_selfloss(p, overlap_min, overlap_max, p3_dist_max)

    # Step 2: Build global binding landscape (tail2weight)
    tail2weight: dict[str, float] = defaultdict(float)
    for primer in primers:
        for tail, weight in primer2tailrc_weight[primer]:
            tail2weight[tail] += weight

    # Step 3: Calculate interaction loss
    interaction_loss = sum(
        tail2weight[tail] * score
        for primer in primers
        for tail, score in primer2tail_score[primer]
    )

    # Step 4: Calculate self-dimerization loss
    self_loss = sum(primer2selfloss[primer] for primer in primers)

    # Step 5: Return total loss
    return interaction_loss + self_loss


if __name__ == "__main__":
    print("=" * 70)
    print("SADDLE Loss Calculation Examples")
    print("=" * 70)
    print()

    # Example 1: Basic usage with a small primer set
    print("Example 1: Basic usage")
    print("-" * 70)
    primers = [
        "ATGCGATCGATCG",
        "CGATCGATCGCAT",
        "GCTAGCTAGCTAG",
    ]
    total_loss = saddle_loss_naive(primers)
    print(f"Primers: {len(primers)} sequences")
    print(f"Total SADDLE loss: {total_loss:,.2f}")
    print()

    # Example 2: Comparing different primer sets
    print("Example 2: Comparing primer sets")
    print("-" * 70)
    set_a = ["ATGCGATCGATCG", "CGATCGATCGCAT"]
    set_b = ["AAAAAAAAAAAAA", "TTTTTTTTTTTT"]
    set_c = ["GCGCGCGCGCGC", "CGCGCGCGCGCG"]

    loss_a = saddle_loss_naive(set_a)
    loss_b = saddle_loss_naive(set_b)
    loss_c = saddle_loss_naive(set_c)

    print(f"Set A (mixed sequence):      {loss_a:>15,.2f}")
    print(f"Set B (homopolymer A/T):     {loss_b:>15,.2f}")
    print(f"Set C (alternating GC):      {loss_c:>15,.2f}")
    print()
    best = min([("A", loss_a), ("B", loss_b), ("C", loss_c)], key=lambda x: x[1])
    print(f"Best primer set: {best[0]} (lower loss is better)")
    print()

    # Example 3: Effect of parameter tuning
    print("Example 3: Parameter sensitivity")
    print("-" * 70)
    test_primers = ["ATGCGATCGATCG", "CGATCGATCGCAT"]

    loss_default = saddle_loss_naive(test_primers)
    loss_strict = saddle_loss_naive(
        test_primers, overlap_min=6, overlap_max=10, p3_dist_max=1
    )
    loss_lenient = saddle_loss_naive(
        test_primers, overlap_min=3, overlap_max=15, p3_dist_max=3
    )

    print(f"Default  (overlap=4-12, p3_dist=2): {loss_default:>15,.2f}")
    print(f"Strict   (overlap=6-10, p3_dist=1): {loss_strict:>15,.2f}")
    print(f"Lenient  (overlap=3-15, p3_dist=3): {loss_lenient:>15,.2f}")
    print()
    print("Note: Stricter parameters focus on stronger, more problematic")
    print("      interactions, while lenient parameters capture weaker ones.")
    print()

    # Example 4: Scaling with primer pool size
    print("Example 4: Scaling with pool size")
    print("-" * 70)
    import random
    random.seed(42)

    def random_primer(length: int = 20) -> str:
        """Generate a random DNA primer."""
        return "".join(random.choice("ACGT") for _ in range(length))

    for pool_size in [2, 5, 10, 20]:
        pool = [random_primer() for _ in range(pool_size)]
        loss = saddle_loss_naive(pool)
        avg_loss = loss / pool_size if pool_size > 0 else 0
        print(f"Pool size {pool_size:>2}: Total={loss:>15,.2f}  Avg/primer={avg_loss:>12,.2f}")

    print()
    print("=" * 70)
    print("Examples complete!")
    print("=" * 70)
