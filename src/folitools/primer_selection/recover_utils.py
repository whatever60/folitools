"""
Species-agnostic gene list simplifier with principled “pseudo_like” filtering and cautious family handling.

Core ideas
----------
1) Name-based pseudo_like vs non-pseudo_like
   - Pseudo-like strings are considered uninformative by *name* alone (placeholders, RIKEN clones, pseudogenes, etc.).
   - Non-pseudo_like are potentially informative and should be preserved.
   - If ≥1 non-pseudo_like symbol exists in the list, **drop all pseudo_like** symbols.
   - If all are pseudo_like, keep them all (order-stable, de-duplicated).

   Pseudo_like patterns (examples):
     • Ensembl IDs (ENSG..., ENSMUSG...), mouse Gm####, RIK (…Rik)
     • lnc-like placeholders (OS/AS/LINC/AC/AL/RP clones)
     • explicit pseudo/related suffixes (-ps, -rs, RNAxxSP#)
     • digits+P at the end (e.g., DEFA10P, CLCA3P)
     • <digit>P<digits> (e.g., KRT8P1), with a safety exception list to avoid real genes like FOXP1, SPP1
     • paralog-like tails where base ends with a digit: <base><lowercase><digits> (e.g., Klk1b1)

2) Readthrough symbols (A-B)
   - A-B is uninformative only when A or B also appears in the same list.
   - Otherwise, keep A-B as informative.

3) Family handling (no over-collapsing)
   - If a **base** token exists, it “eats” its trivial-tail paralogs: collapse to the base (e.g., Klk1 eats Klk1b1, Klk1b3).
   - If the **base is absent**, **do not** aggregate/merge family members (e.g., Selenbp1/Selenbp2 stay separate;
     KRT8P1/KRT8P2 also stay separate). This avoids silently picking an arbitrary representative when the parent is missing.

Rationale for removing earlier ad-hoc strategies
------------------------------------------------
- The prior DEFA-specific collapse and toggleable “aggressive …P<digits>” mode were removed to reduce surprises and
  favor a single transparent rule set. Those strategies either overfit to specific loci (DEFA) or risked collapsing
  bona fide genes (e.g., FOXP1) by name alone. The current approach yields predictable behavior across species and loci.
"""

from collections import OrderedDict
from dataclasses import dataclass
import re


# ---------- compiled patterns (symbol cleanup / pseudo_like) ----------
# ENSG..., ENSMUSG..., optional version
_ENS_GENE_RE = re.compile(r"ENS[A-Z]*G\d{11}(?:\.\d+)?$")
# mouse gene models: Gm10358, Gm3839, ...
_GENE_MODEL_RE = re.compile(r"Gm\d+(?:[A-Za-z0-9]+)?$")
# RIKEN cDNA placeholders (mouse): 1700003D09Rik
_RIK_RE = re.compile(r".*Rik$")
_LNC_LIKE_PATTERNS = (
    re.compile(r".*OS\d*$"),  # MMP24OS, ...OS2
    re.compile(r".*-AS\d*$"),  # FOXP3-AS1, ...-AS2
    re.compile(r".*-AS$"),  # ...-AS
    re.compile(r"^LINC\d+[A-Z]*$"),  # LINC01234
    re.compile(r"^(AC|AL)\d{3,}(?:\.\d+)?$"),  # AC123456.1 / AL123456.1
    re.compile(r"^RP\d+-\d+(?:\.\d+)?$"),  # RP11-...
)
_READTHROUGH_RE = re.compile(r"([A-Za-z0-9]+)-([A-Za-z0-9]+)$")

# Pseudogene/related-sequence patterns
_PSEUDOGENE_DASH_PS_RE = re.compile(r".*-ps\d*$", re.IGNORECASE)  # Clca4c-ps, Defa-ps16
_RELATED_SEQ_DASH_RS_RE = re.compile(r".*-rs\d*$", re.IGNORECASE)  # Rn18s-rs5
_RNA_RRNA_PSEUDO_RE = re.compile(r"RNA\d+SP\d+$")  # RNA18SP5 (human rRNA pseudogene)
_DIGITS_THEN_P_END = re.compile(r".*\d+P$")  # DEFA10P, CLCA3P
# Tighten to require a digit before 'P' to avoid false positives (e.g., MMP24):
_P_THEN_DIGITS_END = re.compile(r".*\dP\d+$")  # ...<digit>P<digits> (e.g., KRT8P1)
_PGENE_PTAIL_RE = re.compile(r"^(?P<base>.+\d)P(?P<num>\d+)$")  # base ends with a digit

# Exceptions to avoid false positives for "...P<digits>" (real genes)
_P_SUFFIX_EXCEPTIONS = {
    "FOXP1",
    "FOXP2",
    "FOXP3",
    "FOXP4",
    "GBP1",
    "ZBP1",
    "GTPBP1",
    "SPP1",
    "DUSP1",
    "AKAP1",
    "TRAP1",
    "RAMP1",
}

# mitochondrial gene prefix (informative)
_MT_PREFIX_RE = re.compile(r"^(?:MT|mt)-")

# ---------- compiled patterns (family collapsing) ----------
# trivial tail for paralog-like variants (lowercase+digits OR digits+lowercase) with optional '-ps'
_FAMILY_TAIL_RE = re.compile(
    r"(?:(?P<lp>[a-z]{1,3})(?P<num>\d+)|(?P<num2>\d+)(?P<ls>[a-z]{1,2}))(?P<ps>-ps)?$"
)


@dataclass(frozen=True)
class TokenDecision:
    """Decision record for a single input token."""

    token: str
    kept: bool
    reason: str


def _is_ensembl_gene_id(symbol: str) -> bool:
    """Return True if symbol looks like an Ensembl gene ID (human/mouse), else False."""
    return bool(_ENS_GENE_RE.fullmatch(symbol))


def _is_gene_model(symbol: str) -> bool:
    """Return True for mouse Gm#### gene models (placeholder ncRNA/protein-coding predictions)."""
    return bool(_GENE_MODEL_RE.fullmatch(symbol))


def _readthrough_parts(symbol: str) -> tuple[str, str] | None:
    """Return (left, right) if `symbol` looks like a simple readthrough A-B; otherwise None."""
    m = _READTHROUGH_RE.fullmatch(symbol)
    if not m:
        return None
    left, right = m.group(1), m.group(2)
    return (left, right) if (left and right) else None


def _is_paralog_like_pseudo(symbol: str) -> bool:
    """Heuristic: symbols like Klk1b1 are pseudo_like when base ends with a digit and tail is [a-z]+digits."""
    m = _FAMILY_TAIL_RE.search(symbol)
    if not m:
        return False
    tail = m.group(0)
    base = symbol[: -len(tail)]
    lp = m.group("lp") or ""
    num = m.group("num") or ""
    # Require lowercase-letter prefix + digits (not digits+lowercase), and base ending with a digit
    return bool(lp and num and base and base[-1].isdigit())


def _is_pseudo_like(symbol: str) -> bool:
    """Name-based heuristic for 'pseudo-like' (uninformative) tokens.

    Args:
        symbol: A gene-like token.

    Returns:
        True if the symbol is considered pseudo_like (uninformative) by naming; else False.
    """
    if _MT_PREFIX_RE.match(symbol):
        return False  # mitochondrial genes are informative
    if _is_ensembl_gene_id(symbol) or _is_gene_model(symbol):
        return True
    if _RIK_RE.fullmatch(symbol):
        return True
    if any(p.fullmatch(symbol) for p in _LNC_LIKE_PATTERNS):
        return True
    if _PSEUDOGENE_DASH_PS_RE.fullmatch(symbol) or _RELATED_SEQ_DASH_RS_RE.fullmatch(
        symbol
    ):
        return True
    if _RNA_RRNA_PSEUDO_RE.fullmatch(symbol):
        return True
    if _DIGITS_THEN_P_END.fullmatch(symbol):
        return True
    if _P_THEN_DIGITS_END.fullmatch(symbol) and symbol not in _P_SUFFIX_EXCEPTIONS:
        return True
    # Paralog-like (e.g., Klk1b1), treated as pseudo when base ends in a digit
    if _is_paralog_like_pseudo(symbol):
        return True
    return False


def _token_decisions(tokens: list[str]) -> list[TokenDecision]:
    """Classify tokens as pseudo_like vs non-pseudo_like and handle readthrough context.

    Args:
        tokens: Order-stable, de-duplicated symbols to classify.

    Returns:
        A list of TokenDecision with per-token keep/drop and reason.
    """
    token_set = set(tokens)
    decisions: list[TokenDecision] = []

    # First, classify pseudo_like by name
    pseudo_flags = {t: _is_pseudo_like(t) for t in tokens}
    any_non_pseudo = any(not pseudo_flags[t] for t in tokens)

    for tok in tokens:
        # Always informative: mitochondrial prefix
        if _MT_PREFIX_RE.match(tok):
            decisions.append(
                TokenDecision(tok, kept=True, reason="informative (mitochondrial)")
            )
            continue

        # Readthrough handling (still informative if components absent)
        parts = _readthrough_parts(tok)
        if parts is not None:
            left, right = parts
            if left in token_set or right in token_set:
                decisions.append(
                    TokenDecision(
                        tok, kept=False, reason="dropped: readthrough_component_present"
                    )
                )
            else:
                decisions.append(
                    TokenDecision(
                        tok,
                        kept=True,
                        reason="informative (readthrough_no_component_present)",
                    )
                )
            continue

        # Pseudo_like dropping rule
        if any_non_pseudo and pseudo_flags[tok]:
            decisions.append(
                TokenDecision(
                    tok,
                    kept=False,
                    reason="dropped: pseudo_like_with_real_symbols_present",
                )
            )
            continue

        # Default: keep (either non-pseudo_like, or all are pseudo_like)
        decisions.append(TokenDecision(tok, kept=True, reason="kept"))

    # If nothing was kept (unlikely), keep all
    if not any(d.kept for d in decisions):
        return [
            TokenDecision(d.token, kept=True, reason="kept: fallback_all_pseudo")
            for d in decisions
        ]
    return decisions


# ---------- family collapsing ----------
@dataclass(frozen=True)
class _TailParts:
    """Parsed trivial tail components for family collapsing."""

    letter_prefix: str
    numeric: int
    letter_suffix: str
    has_ps: bool


def _parse_paralog_tail(symbol: str) -> tuple[str, _TailParts] | None:
    """Parse paralog-like tails with lowercase letters and digits; base must be ≥4 chars.

    Args:
        symbol: Gene symbol to test.

    Returns:
        (base, tail_parts) if a trivial tail is found; else None.
    """
    m = _FAMILY_TAIL_RE.search(symbol)
    if not m:
        return None
    tail = m.group(0)
    base = symbol[: -len(tail)]
    if len(base) < 4:
        return None
    lp = m.group("lp") or ""
    ls = m.group("ls") or ""
    num_str = m.group("num") or m.group("num2")
    num = int(num_str) if num_str is not None else 0
    has_ps = bool(m.group("ps"))
    return base, _TailParts(
        letter_prefix=lp, numeric=num, letter_suffix=ls, has_ps=has_ps
    )


def _parse_pseudogene_ptail(symbol: str) -> tuple[str, int] | None:
    """Parse pseudogene families of the form BASE + 'P' + digits (e.g., KRT8P1).

    Args:
        symbol: Gene symbol to test.

    Returns:
        (base, number) if matches; else None.
    """
    m = _PGENE_PTAIL_RE.match(symbol)
    if not m:
        return None
    base = m.group("base")
    if len(base) < 4:
        return None
    return base, int(m.group("num"))


def _family_sort_key(tp: _TailParts, original_index: int) -> tuple:
    """Sort key for choosing a representative within a family (only used when base exists).

    Args:
        tp: Tail parts (letters, numeric, etc.).
        original_index: Position in the original kept list.

    Returns:
        A tuple used to sort family members deterministically.
    """
    return (tp.letter_prefix, tp.numeric, tp.letter_suffix, tp.has_ps, original_index)


def _collapse_families(kept_tokens: list[str]) -> tuple[list[str], dict[str, str]]:
    """Collapse families with trivial suffix variants (base-aware, no base → no collapse).

    - If base token exists: collapse *all* family members to the base.
    - If base token does NOT exist: **do not** aggregate family members (leave as-is).

    Args:
        kept_tokens: Output of the stage-1 cleanup, order-stable and de-duplicated.

    Returns:
        (new_kept_tokens, dropped_reason_map)
    """
    if not kept_tokens:
        return kept_tokens, {}

    token_set = set(kept_tokens)

    by_base_paralog: dict[str, list[tuple[int, str, _TailParts]]] = {}
    by_base_pgene: dict[str, list[tuple[int, str, int]]] = {}

    for idx, tok in enumerate(kept_tokens):
        parsed_para = _parse_paralog_tail(tok)
        if parsed_para:
            base, tail_parts = parsed_para
            by_base_paralog.setdefault(base, []).append((idx, tok, tail_parts))
        parsed_pgene = _parse_pseudogene_ptail(tok)
        if parsed_pgene:
            base, num = parsed_pgene
            by_base_pgene.setdefault(base, []).append((idx, tok, num))

    to_drop: dict[str, str] = {}

    # Paralog families
    for base, members in by_base_paralog.items():
        if base in token_set:
            for _, tok, _ in members:
                to_drop[tok] = f"collapsed_into_base: base={base}"
        # base absent → no aggregation

    # Pseudogene P<digits> families
    for base, members in by_base_pgene.items():
        if base in token_set:
            for _, tok, _ in members:
                to_drop[tok] = f"collapsed_into_base: base={base}"
        # base absent → no aggregation

    new_kept = [t for t in kept_tokens if t not in to_drop]
    return new_kept, to_drop


# ---------- public API ----------
def _simplify_gene_list_internal(
    gene_names: list[str],
    *,
    with_audit: bool = False,
    collapse_families: bool = True,
) -> tuple[list[str], dict[str, str]]:
    """Internal implementation for gene list simplification per 'pseudo_like' rules.

    Args:
        gene_names: List of raw gene symbols/IDs.
        with_audit: If True, include dropped reasons in return value.
        collapse_families: If True, apply the family collapsing pass.

    Returns:
        (kept_list, dropped_reason_map). When with_audit=False, dropped_reason_map is {}.
    """
    tokens: list[str] = list(
        OrderedDict((t.strip(), None) for t in gene_names if t.strip()).keys()
    )

    decisions = _token_decisions(tokens)
    kept_tokens = [d.token for d in decisions if d.kept]
    dropped: dict[str, str] = (
        {d.token: d.reason for d in decisions if not d.kept} if with_audit else {}
    )

    if collapse_families and kept_tokens:
        kept_tokens, fam_dropped = _collapse_families(kept_tokens)
        if with_audit and fam_dropped:
            dropped.update(fam_dropped)

    return kept_tokens, dropped


def simplify_gene_list(
    gene_names: list[str],
    *,
    collapse_families: bool = True,
) -> list[str]:
    """Simplify a list of gene symbols/IDs per 'pseudo_like' rules and base-aware family collapsing.

    Args:
        gene_names: List of gene symbols/IDs.
        collapse_families: If True, apply the family collapsing pass (base required to collapse).

    Returns:
        List of kept tokens in first-seen order (post-collapsing).
    """
    kept, _ = _simplify_gene_list_internal(
        gene_names, with_audit=False, collapse_families=collapse_families
    )
    return kept


def simplify_gene_list_with_audit(
    gene_names: list[str],
    *,
    collapse_families: bool = True,
) -> tuple[list[str], dict[str, str]]:
    """Simplify a list of gene symbols/IDs and also return a map of dropped tokens → reason.

    Args:
        gene_names: List of gene symbols/IDs.
        collapse_families: If True, apply the family collapsing pass (base required to collapse).

    Returns:
        (kept_list, dropped_reason_map).
    """
    return _simplify_gene_list_internal(
        gene_names, with_audit=True, collapse_families=collapse_families
    )


# ---------- tests / examples ----------
if __name__ == "__main__":

    def expect(
        inp: list[str], out: list[str], *, collapse_families: bool = True
    ) -> None:
        got = simplify_gene_list(inp, collapse_families=collapse_families)
        assert got == out, f"{inp} → {got}, expected {out}"

    # Stage-1 sanity
    # A readthrough is dropped if a component is present
    expect(["KLRC4-KLRK1", "KLRK1"], ["KLRK1"])
    # A readthrough is kept if alone
    expect(["IFNAR2-IL10RB"], ["IFNAR2-IL10RB"])

    # All pseudo_like -> keep all (order-stable)
    expect(["ENSG00000173366", "Gm3839"], ["ENSG00000173366", "Gm3839"])

    # RIK placeholders dropped when a real symbol is present
    expect(["Ppp1r1b", "1700003D09Rik"], ["Ppp1r1b"])

    # lnc-like / models
    expect(["Gm10358", "Gapdh", "Gm3839"], ["Gapdh"])
    expect(["MMP24", "MMP24OS"], ["MMP24"])
    expect(["FOXP3-AS1", "FOXP3"], ["FOXP3"])

    # Mitochondrial
    # ENSG0000 is not a real Ensembl id, so both are kept.
    expect(["MT-CO3", "ENSG0000"], ["MT-CO3", "ENSG0000"])

    # Pseudogenes vs real
    # DEFA10P is pseudo_like relative to DEFA1
    expect(["DEFA10P", "DEFA1"], ["DEFA1"])
    # CLCA3P is pseudo_like relative to CLCA3
    expect(["CLCA3P", "CLCA3"], ["CLCA3"])
    # RNA18SP5 is pseudo_like relative to RNA18S
    expect(["RNA18SP5", "RNA18S"], ["RNA18S"])
    # HMGB1P1 is pseudo_like relative to HMGB1
    expect(["HMGB1P1", "HMGB1"], ["HMGB1"])

    # Non-pseudo paralogs WITHOUT base present: keep both (Selenbp family)
    expect(["Selenbp1", "Selenbp2"], ["Selenbp1", "Selenbp2"])

    # Pseudogene-only families WITHOUT base: DO NOT aggregate — keep all
    expect(["KRT8P1", "KRT8P2"], ["KRT8P1", "KRT8P2"])

    # Mixed families + outsiders: real symbol drops pseudo_like
    expect(["Klk1b1", "Klk1b2", "Klk2"], ["Klk2"])

    # Base present → collapse to base
    expect(
        [
            "Klk1",
            "Klk1b1",
            "Klk1b11",
            "Klk1b14-ps",
            "Klk1b21",
            "Klk1b22",
            "Klk1b3",
            "Klk1b4",
            "Klk1b5",
            "Klk1b8",
            "Klk1b9",
        ],
        ["Klk1"],
    )

    # Do not mix with nontrivial tails (uppercase letters), and base absent → keep all
    expect(["H2BC12", "H2BC12L", "H2BC13"], ["H2BC12", "H2BC12L", "H2BC13"])

    print("All tests passed ✓")
