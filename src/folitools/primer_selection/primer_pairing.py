"""
primer_pairing.py

Fast complementary subsequence extraction for DNA primers (A/C/G/T only).

This module provides two main functions:

1) pair_primers(...)
   - Given two primers in 5'->3' orientation, find complementary subsequences
     of length in [len_min, len_max], where:
         subseq2 == reverse_complement(subseq1)
   - Returns distances to each primer's 3' end (0 means the subseq ends at 3').

2) pair_primer_pool(...)
   - Given a pool of primers, stream all complementary hits between primers
     without explicit O(M^2) pair scanning.
   - Builds a global inverted index once, then joins each primer's reverse-
     complement k-mers against that index.

Key implementation ideas:
- 2-bit encoding of DNA bases: A=0, C=1, G=2, T=3
- Rolling k-mer integers: update window value in O(1) as you slide by 1 base
- Precompute reverse-complement once per primer (idea 2)

Important constraints:
- Only A/C/G/T are allowed. N is rejected.
- Primers are treated as given in 5'->3' orientation.
"""

from collections import defaultdict
from collections.abc import Iterator
from dataclasses import dataclass
import tempfile
from pathlib import Path

import polars as pl
from tqdm.auto import tqdm


_BASE_TO_CODE: dict[str, int] = {"A": 0, "C": 1, "G": 2, "T": 3}
_CODE_TO_BASE: str = "ACGT"


@dataclass(frozen=True, slots=True)
class _Occurrence:
    """
    A single k-mer occurrence inside one primer.

    Attributes:
        primer_id: Index of the primer in the input primer list.
        start_5p: 0-based start position of the k-mer in the primer (5'->3' string).
        dist_to_3p: Distance from the k-mer end to the primer's 3' end.
                   0 means the k-mer ends at the 3' base.
    """

    primer_id: int
    start_5p: int
    dist_to_3p: int


def _validate_len_bounds(len_min: int, len_max: int) -> None:
    if len_min <= 0 or len_max <= 0 or len_min > len_max:
        raise ValueError(f"Invalid len_min/len_max: {len_min=}, {len_max=}")


def _validate_max_dist_to_3p(max_dist_to_3p: int | None) -> None:
    if max_dist_to_3p is None:
        return
    if max_dist_to_3p < 0:
        raise ValueError(
            f"max_dist_to_3p must be >= 0 (or None), got {max_dist_to_3p!r}"
        )


def _encode_acgt_to_bytes(seq: str, *, label: str) -> bytes:
    """
    Encode a DNA string (A/C/G/T only) into bytes of integers 0..3.

    Args:
        seq: DNA string.
        label: Used in error messages.

    Returns:
        bytes of length len(seq), with each byte in {0,1,2,3}.
    """
    s = seq.upper()
    bad = {b for b in s if b not in _BASE_TO_CODE}
    if bad:
        raise ValueError(f"{label} has invalid bases {sorted(bad)}: {seq!r}")
    return bytes(_BASE_TO_CODE[b] for b in s)


def _reverse_complement_codes_bytes(codes: bytes) -> bytes:
    """
    Reverse-complement in 2-bit code space.

    With encoding A=0,C=1,G=2,T=3, complement is code ^ 3.
    """
    return bytes((c ^ 3) for c in reversed(codes))


def reverse_complement(dna_5to3: str) -> str:
    """
    Reverse-complement a DNA string (A/C/G/T only).

    Args:
        dna_5to3: DNA string.

    Returns:
        Reverse-complement string in 5'->3' orientation.
    """
    codes = _encode_acgt_to_bytes(dna_5to3, label="dna_5to3")
    rc_codes = _reverse_complement_codes_bytes(codes)
    return "".join(_CODE_TO_BASE[c] for c in rc_codes)


def _iter_rolling_kmer_ints(
    codes: bytes,
    k: int,
) -> Iterator[tuple[int, int]]:
    """
    Yield (start, kmer_int) for every k-mer window in a code sequence.

    Uses rolling update:
        val_next = ((val << 2) | new_code) & mask

    Args:
        codes: bytes of 0..3.
        k: k-mer length.

    Yields:
        (start_index, kmer_int) for start in [0, n-k].
    """
    n = len(codes)
    if k <= 0 or k > n:
        return

    mask = (1 << (2 * k)) - 1
    val = 0
    for t in range(k):
        val = (val << 2) | codes[t]

    yield 0, val
    for start in range(1, n - k + 1):
        val = ((val << 2) | codes[start + k - 1]) & mask
        yield start, val


def _start_min_for_3p_window(n: int, k: int, max_dist_to_3p: int | None) -> int:
    """
    Compute the minimum 0-based start position such that dist_to_3p <= W.

    dist_to_3p = n - (start + k) <= W
        start >= n - k - W

    Args:
        n: sequence length.
        k: window length.
        max_dist_to_3p: W, or None to disable restriction.

    Returns:
        start_min (clamped to >= 0).
    """
    if max_dist_to_3p is None:
        return 0
    return max(0, n - k - max_dist_to_3p)


def pair_primers(
    seq1_5to3: str,
    seq2_5to3: str,
    len_min: int = 4,
    len_max: int = 8,
    max_dist_to_3p: int | None = None,
) -> list[tuple[str, str, int, int]]:
    """
    Find complementary subsequences between two primers (5'->3').

    A hit is defined as:
        subseq2 == reverse_complement(subseq1)

    Distances are to the 3' end:
        dist_to_3p = number of bases between subseq end and primer 3' end
        so subseq ending exactly at the 3' base has dist_to_3p = 0.

    Optional 3' window restriction:
        If max_dist_to_3p = W is provided, only consider k-mers whose dist_to_3p <= W
        on BOTH primers.

    Args:
        seq1_5to3: Primer 1 in 5'->3' (A/C/G/T only).
        seq2_5to3: Primer 2 in 5'->3' (A/C/G/T only).
        len_min: Minimum subseq length (default 4).
        len_max: Maximum subseq length (default 8).
        max_dist_to_3p: If not None, restrict to windows ending within W bases of 3'
                        (default None: no restriction).

    Returns:
        List of tuples:
            (subseq_in_seq1, subseq_in_seq2, dist1_to_3p, dist2_to_3p)

        If you need 0-based start positions:
            start = len(seq) - dist_to_3p - k
    """
    _validate_len_bounds(len_min, len_max)
    _validate_max_dist_to_3p(max_dist_to_3p)

    s1 = seq1_5to3.upper()
    s2 = seq2_5to3.upper()
    n1 = len(s1)
    n2 = len(s2)

    codes1 = _encode_acgt_to_bytes(s1, label="seq1_5to3")
    codes2 = _encode_acgt_to_bytes(s2, label="seq2_5to3")
    codes1_rc = _reverse_complement_codes_bytes(codes1)

    results: list[tuple[str, str, int, int]] = []
    max_k = min(len_max, n1, n2)

    for k in range(len_min, max_k + 1):
        # Index k-mers in seq2 by integer value -> list of starts (restricted to 3' window if requested).
        start2_min = _start_min_for_3p_window(n2, k, max_dist_to_3p)
        kmer2_to_starts: dict[int, list[int]] = defaultdict(list)

        for start2, v2 in _iter_rolling_kmer_ints(codes2, k):
            if start2 < start2_min:
                continue
            kmer2_to_starts[v2].append(start2)

        # Query seq1 windows (restricted to 3' window if requested):
        start1_min = _start_min_for_3p_window(n1, k, max_dist_to_3p)

        # We fetch RC(subseq1) by using the rolling k-mer ints on rc(seq1) codes.
        # For subseq1 starting at start1 in seq1, RC(subseq1) starts at:
        #     p = n1 - (start1 + k)
        # and codes1_rc[p:p+k] is exactly RC(subseq1).
        # Rather than building a full list, we just iterate all starts on codes1_rc
        # and look up the relevant p each time; since n1 is small (~25), this is fast.
        # To avoid re-scanning codes1_rc for each start1, we precompute all (start, val)
        # into a list for the current k (still tiny memory for primers).
        rc_kmers: list[int] = [0] * (n1 - k + 1)
        for p, v in _iter_rolling_kmer_ints(codes1_rc, k):
            rc_kmers[p] = v

        for start1 in range(start1_min, n1 - k + 1):
            dist1 = n1 - (start1 + k)
            p = dist1  # same as n1 - (start1 + k)
            needed = rc_kmers[p]
            starts2 = kmer2_to_starts.get(needed)
            if not starts2:
                continue

            sub1 = s1[start1 : start1 + k]
            for start2 in starts2:
                dist2 = n2 - (start2 + k)
                sub2 = s2[start2 : start2 + k]
                results.append((sub1, sub2, dist1, dist2))

    return results


def pair_primer_pool(
    primers_5to3: list[str],
    len_min: int = 4,
    len_max: int = 8,
    include_self: bool = True,
    max_dist_to_3p: int | None = None,
) -> Iterator[tuple[int, int, str, str, int, int]]:
    """
    Stream complementary hits between primers drawn from one pool.

    Yield format:
        (primer1_id, primer2_id, subseq1, subseq2, dist1_to_3p, dist2_to_3p)

    Optional 3' window restriction:
        If max_dist_to_3p = W is provided, only consider k-mers whose dist_to_3p <= W
        BOTH when building the pool index and when querying each primer.

    This restriction usually:
      - dramatically reduces index size (memory)
      - dramatically reduces hit counts for small k (time/output)
      - matches biology (3' complementarity is most extension-relevant)

    Args:
        primers_5to3: List of primers (A/C/G/T only), each in 5'->3' orientation.
        len_min: Minimum subseq length (default 4).
        len_max: Maximum subseq length (default 8).
        include_self: If False, do not yield hits where primer1_id == primer2_id.
        max_dist_to_3p: If not None, restrict windows ending within W bases of 3'
                        (default None: no restriction).

    Yields:
        (primer1_id, primer2_id, subseq1, subseq2, dist1_to_3p, dist2_to_3p)
    """
    _validate_len_bounds(len_min, len_max)
    _validate_max_dist_to_3p(max_dist_to_3p)

    primers_u = [p.upper() for p in primers_5to3]
    if not primers_u:
        return
        yield  # pragma: no cover

    primers_codes = [
        _encode_acgt_to_bytes(p, label=f"primers_5to3[{i}]")
        for i, p in enumerate(primers_u)
    ]
    lens = [len(p) for p in primers_u]
    global_max_len = min(len_max, min(lens))

    # Build inverted index once:
    #   index_by_k[k][kmer_int] -> list of occurrences across the pool.
    index_by_k: dict[int, dict[int, list[_Occurrence]]] = {}

    for k in range(len_min, global_max_len + 1):
        idx: dict[int, list[_Occurrence]] = defaultdict(list)

        for pid, codes in enumerate(primers_codes):
            n = len(codes)
            if n < k:
                continue

            start_min = _start_min_for_3p_window(n, k, max_dist_to_3p)
            for start, kmer_int in _iter_rolling_kmer_ints(codes, k):
                if start < start_min:
                    continue
                dist_to_3p = n - (start + k)
                idx[kmer_int].append(_Occurrence(pid, start, dist_to_3p))

        index_by_k[k] = idx

    # Query each primer against the global index.
    for pid1, (s1, codes1) in enumerate(
        zip(tqdm(primers_u), primers_codes, strict=True)
    ):
        n1 = len(s1)
        rc_codes1 = _reverse_complement_codes_bytes(codes1)

        for k in range(len_min, min(global_max_len, n1) + 1):
            idx = index_by_k[k]

            # Precompute rolling ints for rc(seq1) for this k into a small list.
            rc_kmers: list[int] = [0] * (n1 - k + 1)
            for p, v in _iter_rolling_kmer_ints(rc_codes1, k):
                rc_kmers[p] = v

            start1_min = _start_min_for_3p_window(n1, k, max_dist_to_3p)

            for start1 in range(start1_min, n1 - k + 1):
                dist1 = n1 - (start1 + k)
                p = dist1
                key = rc_kmers[p]

                # Example (0-based indexing):
                #
                # s1     = A C G T A C C A
                # idx      0 1 2 3 4 5 6 7
                #              ^
                #              start1 = 2, k = 3  =>  s1[start1:start1+k] = "GTA"
                # rc(s1) = T G G T A C G T
                # idx      0 1 2 3 4 5 6 7
                #                ^
                #                p = n1 - (start1 + k) = 8 - (2 + 3) = 3
                #                rc(s1)[p:p+k] = rc(s1)[3:6] = "TAC" == RC("GTA")
                # So `key = kmers_rc[p]` is the rolling-kmer key representing the window
                # rc(s1)[p:p+k] (i.e., the reverse-complement of s1[start1:start1+k]).
                # RC(subseq1) maps to a window in rc(seq1) starting at p.

                hits = idx.get(key)
                if not hits:
                    continue

                sub1 = s1[start1 : start1 + k]
                for occ in hits:
                    if (not include_self) and (occ.primer_id == pid1):
                        continue
                    sub2 = primers_u[occ.primer_id][occ.start_5p : occ.start_5p + k]
                    yield (pid1, occ.primer_id, sub1, sub2, dist1, occ.dist_to_3p)


def collect_pair_primer_pool_results(iterator: object) -> pl.DataFrame:
    """
    Collect yielded pairwise primer-pool results into a Polars DataFrame using batched Parquet writes.

    This is memory-friendly for large iterators: it writes records in batches to Parquet part-files,
    then reads them back via `scan_parquet(...).collect()`.

    Args:
        iterator: An iterator that yields 6-tuples in the form:
            (primer_id_1, primer_id_2, subseq_1, subseq_2, dist_to_3p_1, dist_to_3p_2)

    Returns:
        A Polars DataFrame with columns:
            - primer_id_1 (UInt32)
            - primer_id_2 (UInt32)
            - subseq_1 (Utf8)
            - subseq_2 (Utf8)
            - dist_to_3p_1 (UInt32)
            - dist_to_3p_2 (UInt32)
    """
    schema = [
        ("primer_id_1", pl.UInt32),
        ("primer_id_2", pl.UInt32),
        ("subseq_1", pl.Utf8),
        # ("subseq_2", pl.Utf8),
        ("dist_to_3p_1", pl.UInt32),
        ("dist_to_3p_2", pl.UInt32),
    ]

    batch_size = 200_000
    batch: list[tuple[int, int, str, int, int]] = []
    part = 0
    wrote_any = False

    def _empty_df() -> pl.DataFrame:
        return pl.DataFrame(
            {name: pl.Series(name, [], dtype=dtype) for name, dtype in schema}
        )

    with tempfile.TemporaryDirectory(prefix="primer_pool_hits_") as tmp:
        out_dir = Path(tmp)

        for rec in iterator:
            try:
                pid1, pid2, sub1, sub2, dist1, dist2 = rec
            except Exception as e:
                raise ValueError(
                    "Expected each yielded record to be a 6-tuple: "
                    "(primer_id_1, primer_id_2, subseq_1, subseq_2, dist_to_3p_1, dist_to_3p_2)"
                ) from e

            # batch.append((pid1, pid2, sub1, sub2, dist1, dist2))
            batch.append((pid1, pid2, sub1, dist1, dist2))

            if len(batch) >= batch_size:
                pl.DataFrame(batch, schema=schema, orient="row").write_parquet(
                    out_dir / f"part-{part:05d}.parquet"
                )
                wrote_any = True
                batch.clear()
                part += 1

        if batch:
            pl.DataFrame(batch, schema=schema, orient="row").write_parquet(
                out_dir / f"part-{part:05d}.parquet"
            )
            wrote_any = True

        if not wrote_any:
            return _empty_df()

        df = pl.scan_parquet(str(out_dir / "part-*.parquet"))

        return (
            df.with_columns(
                pl.col("subseq_1").str.count_matches(r"[CG]").alias("gc_count"),
                pl.col("subseq_1").str.len_chars().alias("len"),
            )
            .with_columns(
                (
                    pl.lit(2) ** (pl.col("len") + pl.col("gc_count"))
                    / (
                        (pl.lit(1) + pl.col("dist_to_3p_1"))
                        * (pl.lit(1) + pl.col("dist_to_3p_2"))
                    )
                    # Since non-self interactions are counted twice (A->B and B->A), we
                    # double the score for self-interactions to keep the total sum consistent.
                    * pl.when(pl.col("primer_id_1").eq(pl.col("primer_id_2")))
                    .then(2)
                    .otherwise(1)
                ).alias("score")
            )
            .collect()
        )


if __name__ == "__main__":
    import time

    import numpy as np

    from .saddle_loss import saddle_loss_naive

    # ----------------------------
    # Example 1: pair two primers
    # ----------------------------
    p1 = "ACGTTGCA"
    p2 = "TGCAACGT"  # reverse-complement of p1
    # p1 = "ATATATAT"
    # p2 = "TATATATA"  # reverse-complement of p1

    hits = pair_primers(p1, p2, len_min=4, len_max=8)
    print("Pairwise hits (showing first 10):")
    for row in hits[:10]:
        sub1, sub2, d1, d2 = row
        k = len(sub1)
        # Recover 0-based start positions from distances if desired:
        start1 = len(p1) - d1 - k
        start2 = len(p2) - d2 - k
        print(
            f"  sub1={sub1} (start={start1}, dist to 3'={d1})  sub2={sub2} "
            f"(start={start2}, dist to 3'={d2})"
        )
    print(f"Total hits: {len(hits)}\n")

    # A quick sanity check on the distance convention:
    # Subsequence ending at the 3' base has dist_to_3p = 0.
    p = "AACCGGTT"
    last4 = p[-4:]
    demo_hits = pair_primers(p, reverse_complement(p), len_min=4, len_max=4)
    # Find the hit that corresponds to the last 4 bases of p (ends at 3').
    for sub1, sub2, d1, d2 in demo_hits:
        if sub1 == last4:
            print("Distance-to-3' sanity check:")
            print(f"  sub1={sub1} ends at 3' => dist1_to_3p={d1} (expected 0)")
            break
    print()

    # ----------------------------
    # Example 2: pair a primer pool
    # ----------------------------
    pool = [
        "ACGTTGCA",
        "TGCAACGT",  # RC of pool[0]
        "AAAACCCC",
        "GGGGTTTT",  # RC of pool[2]
        "ATATATAT",  # repetitive; will create more short-k collisions
    ]

    print("Pool hits (k=6..8), first 20 yielded:")
    gen = pair_primer_pool(pool, len_min=6, len_max=8, include_self=False)
    for i, hit in enumerate(gen):
        if i >= 20:
            break
        pid1, pid2, sub1, sub2, d1, d2 = hit
        print(f"  ({pid1}->{pid2}) {sub1} <-> {sub2} dist to 3' = ({d1}, {d2})")

    # Tip for large pools:
    # - prefer len_min >= 7 or 8
    # - consider restricting to 3' windows if you only care about extension-relevant pairing

    # ----------------------------
    # Random primer generation (fast via numpy if available)
    # ----------------------------
    def make_random_primers(n: int, length: int, seed: int = 123) -> list[str]:
        bases = "ACGT"

        rng = np.random.default_rng(seed)
        arr = rng.integers(0, 4, size=(n, length), dtype=np.uint8)
        # convert each row to a string
        return ["".join(bases[x] for x in row.tolist()) for row in arr]

    primer_len = 25
    n_primers = 10000
    primers_all = make_random_primers(n_primers, primer_len, seed=123)

    t0 = time.perf_counter()
    iterator = pair_primer_pool(primers_all, len_min=4, len_max=8, max_dist_to_3p=2)
    res = collect_pair_primer_pool_results(iterator)
    t1 = time.perf_counter()
    # naive_loss = saddle_loss_naive(
    #     primers_all, overlap_min=4, overlap_max=8, p3_dist_max=2
    # )
    naive_loss = 0
    print(
        f"Timing note: pairing {n_primers} primers (len_min=4) took {t1 - t0:.2f} "
        f"seconds and yielded: {len(res)} hits.\n"
        f"Total score: {res['score'].sum()}\n"
        f"Naive SADDLE loss: {naive_loss}\n"
    )
    print(res.head())
