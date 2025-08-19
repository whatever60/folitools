from collections.abc import Sequence
import random
from typing import TypeVar

T = TypeVar("T")

def choice_except(seq: Sequence[T], i: int, *, rng: random.Random | None = None) -> T:
    """Choose a random element from `seq` excluding the element at index `i`.

    Args:
        seq: The input sequence (e.g., list, tuple).
        i: Index to exclude; negative indices are supported (Python-style).
        rng: Optional random generator to control determinism; defaults to `random`.

    Returns:
        A randomly chosen element from `seq` that is not at position `i`.

    Raises:
        ValueError: If `seq` has fewer than 2 elements, or `i` is out of range.
    """
    n = len(seq)
    if n <= 1:
        raise ValueError("Sequence must have at least two elements to exclude one.")
    j = i if i >= 0 else n + i
    if not 0 <= j < n:
        raise ValueError(f"Index i={i} out of range for sequence of length {n}.")

    # Uniformly sample an index from the n-1 allowable positions without copying.
    randrange = (rng.randrange if rng else random.randrange)
    k = randrange(n - 1)
    if k >= j:
        k += 1
    return seq[k]
