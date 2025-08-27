"""
Species-agnostic gene list simplifier with family-collapsing.

Pipeline (conservative by default):
  1) Symbol cleanup:
     • If ≥1 informative symbols exist, keep ONLY those; otherwise keep ALL tokens (order-stable, de-duplicated).
     • Readthrough 'A-B' is uninformative ONLY when A or B is also present; otherwise it is informative.
     • Pseudogene/related-sequence/gene-model/lnc-like placeholders/plain Ensembl IDs are uninformative.
       - Safe patterns (on by default): mouse '-ps'/'-rs', 'RNA18SP#', '<digits>P' (e.g., DEFA10P),
         'Gm####' gene models, OS/AS/LINC/AC/AL/RP clones, Ensembl gene IDs.
       - Risky human '...P<digits>' (e.g., HMGB1P1) is OFF by default; enable via Rules(aggressive_human_p_suffix=True).

  2) Family collapsing (new):
     • Cluster remaining tokens by a nontrivial base (≥4 chars) + a TRIVIAL tail.
       - TRIVIAL tail := ends with one of:
           (a) [a-z]{1,3} + digits, e.g., 'b1', 'b11', 'ab3'
           (b) digits + [a-z]{0,2}, e.g., '1', '2', '3a'
         optionally followed by '-ps' (e.g., 'b14-ps').
       - Uppercase letters are NOT trivial, so 'H2BC12L' won't collapse with 'H2BC12'.
     • For each family (same base):
         (i)  if the base token itself exists, collapse all family members to the base.
         (ii) else, keep the single earliest member by natural order of the parsed tail
              (letter_prefix, numeric, letter_suffix, non-ps before ps), drop the rest.
     • Families are independent; non-members are left untouched.

I/O:
  - simplify_gene_list(List[str]) -> List[str]
  - simplify_gene_list_with_audit(List[str]) -> (List[str], Dict[str, str])

"""

from collections import OrderedDict
from dataclasses import dataclass
import re


# ---------- compiled patterns (symbol cleanup) ----------
_ENS_GENE_RE = re.compile(
    r"ENS[A-Z]*G\d{11}(?:\.\d+)?$"
)  # ENSG..., ENSMUSG..., optional version
_GENE_MODEL_RE = re.compile(
    r"Gm\d+(?:[A-Za-z0-9]+)?$"
)  # mouse gene models: Gm10358, Gm3839, ...
_LNC_LIKE_PATTERNS = (
    re.compile(r".*OS\d*$"),  # MMP24OS, ...OS2
    re.compile(r".*-AS\d*$"),  # FOXP3-AS1, ...-AS2
    re.compile(r".*-AS$"),  # ...-AS
    re.compile(r"^LINC\d+[A-Z]*$"),  # LINC01234
    re.compile(r"^(AC|AL)\d{3,}(?:\.\d+)?$"),  # AC123456.1 / AL123456.1
    re.compile(r"^RP\d+-\d+(?:\.\d+)?$"),  # RP11-...
)
_READTHROUGH_RE = re.compile(r"([A-Za-z0-9]+)-([A-Za-z0-9]+)$")

# Pseudogene/related-sequence: SAFE patterns valid in both species
_PSEUDOGENE_DASH_PS_RE = re.compile(r".*-ps\d*$", re.IGNORECASE)  # Clca4c-ps, Defa-ps16
_RELATED_SEQ_DASH_RS_RE = re.compile(r".*-rs\d*$", re.IGNORECASE)  # Rn18s-rs5
_RNA_RRNA_PSEUDO_RE = re.compile(r"RNA\d+SP\d+$")  # RNA18SP5 (human rRNA pseudogene)
_DIGITS_THEN_P_END = re.compile(
    r".*\d+P$"
)  # DEFA10P, CLCA3P (digits before terminal 'P')

# Aggressive human-style “…P<digits>” (risky — off by default)
_P_THEN_DIGITS_END = re.compile(
    r".*P\d+$"
)  # HMGB1P1, RPL23AP82, etc. (also matches FOXP1/AKAP1/...!)

# mitochondrial gene prefix (informative)
_MT_PREFIX_RE = re.compile(r"^(?:MT|mt)-")

# ---------- compiled patterns (family collapsing) ----------
# trivial tail: [a-z]{1,3}\d+  OR  \d+[a-z]{0,2}, optional '-ps'
_FAMILY_TAIL_RE = re.compile(
    r"(?:(?P<lp>[a-z]{1,3})(?P<num>\d+)|(?P<num2>\d+)(?P<ls>[a-z]{0,2}))(?P<ps>-ps)?$"
)


@dataclass(frozen=True)
class TokenDecision:
    """Decision record for a single input token."""

    token: str
    kept: bool
    reason: str


@dataclass(frozen=True)
class Rules:
    """Configuration for species-agnostic simplification.

    Args:
        aggressive_human_p_suffix: If True, treat symbols ending with 'P<digits>' as pseudogenes.
            WARNING: This can misclassify bona fide genes (e.g., FOXP1, AKAP1, DUSP1, GBP1, ZBP1, GTPBP1, SPP1, TRAP1).
            Keep False unless you really need to collapse these by symbol alone.
    """

    aggressive_human_p_suffix: bool = False


# ---------- helpers: stage 1 (symbol cleanup) ----------
def _is_ensembl_gene_id(symbol: str) -> bool:
    return bool(_ENS_GENE_RE.fullmatch(symbol))


def _is_gene_model(symbol: str) -> bool:
    return bool(_GENE_MODEL_RE.fullmatch(symbol))


def _is_lnc_like_symbol(symbol: str) -> bool:
    return any(p.fullmatch(symbol) for p in _LNC_LIKE_PATTERNS)


def _readthrough_parts(symbol: str) -> tuple[str, str] | None:
    m = _READTHROUGH_RE.fullmatch(symbol)
    if not m:
        return None
    left, right = m.group(1), m.group(2)
    return (left, right) if (left and right) else None


def _is_pseudogene_or_related(symbol: str, rules: Rules) -> bool:
    """Detect pseudogene/related-sequence by symbol only (species-agnostic, conservative)."""
    if _PSEUDOGENE_DASH_PS_RE.fullmatch(symbol):
        return True
    if _RELATED_SEQ_DASH_RS_RE.fullmatch(symbol):
        return True
    if _RNA_RRNA_PSEUDO_RE.fullmatch(symbol):
        return True
    if _DIGITS_THEN_P_END.fullmatch(symbol):  # e.g., DEFA10P, CLCA3P
        return True
    # Aggressive human pattern (off by default): ...P<digits>
    if rules.aggressive_human_p_suffix and _P_THEN_DIGITS_END.fullmatch(symbol):
        return True
    return False


def _token_decisions(tokens: list[str], rules: Rules) -> list[TokenDecision]:
    """Classify tokens as informative or not; readthroughs are context-aware.

    Args:
        tokens: Order-stable, de-duplicated symbols to classify.
        rules: Behavior toggles for conservative vs. aggressive pseudogene handling.

    Returns:
        A list of TokenDecision with per-token keep/drop and reason.
    """
    token_set = set(tokens)
    decisions: list[TokenDecision] = []

    for tok in tokens:
        # Always informative: mitochondrial prefix (both species)
        if _MT_PREFIX_RE.match(tok):
            decisions.append(
                TokenDecision(tok, kept=True, reason="informative (mitochondrial)")
            )
            continue

        # Hard uninformative classes (symbol-only)
        if _is_ensembl_gene_id(tok):
            decisions.append(
                TokenDecision(tok, kept=False, reason="dropped: ensembl_id")
            )
            continue
        if _is_gene_model(tok):
            decisions.append(
                TokenDecision(tok, kept=False, reason="dropped: gene_model_placeholder")
            )
            continue
        if _is_lnc_like_symbol(tok):
            decisions.append(
                TokenDecision(tok, kept=False, reason="dropped: lnc_like_placeholder")
            )
            continue
        if _is_pseudogene_or_related(tok, rules):
            decisions.append(
                TokenDecision(tok, kept=False, reason="dropped: pseudogene_or_related")
            )
            continue

        # Readthrough is uninformative only if a component appears
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

        # Default informative
        decisions.append(TokenDecision(tok, kept=True, reason="informative"))

    # If nothing was informative, keep everything (order-stable)
    if not any(d.kept for d in decisions):
        return [
            TokenDecision(d.token, kept=True, reason="kept: all_tokens_uninformative")
            for d in decisions
        ]
    return decisions


# ---------- helpers: stage 2 (family collapsing) ----------
@dataclass(frozen=True)
class _TailParts:
    """Parsed trivial tail components for family collapsing."""

    letter_prefix: str
    numeric: int
    letter_suffix: str
    has_ps: bool


def _parse_trivial_tail(symbol: str) -> tuple[str, _TailParts] | None:
    """Split a symbol into (base, tail_parts) if it ends with a trivial tail.

    A nontrivial base must have length ≥ 4.

    Args:
        symbol: Gene symbol to parse.

    Returns:
        (base, _TailParts) if a trivial tail is found and base is long enough; otherwise None.
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


def _family_sort_key(tp: _TailParts, original_index: int) -> tuple:
    """Sort key: (letter_prefix, numeric asc, letter_suffix, non-ps first, original_index)."""
    return (tp.letter_prefix, tp.numeric, tp.letter_suffix, tp.has_ps, original_index)


def _collapse_families(kept_tokens: list[str]) -> tuple[list[str], dict[str, str]]:
    """Collapse families with trivial suffix variants.

    Rules:
      • If base token exists in kept_tokens: collapse all base+tail members to the base.
      • Else: keep the earliest member by natural tail order; drop the rest.
      • Operates independently per base; non-members untouched.

    Args:
        kept_tokens: Output of the stage-1 cleanup, order-stable and de-duplicated.

    Returns:
        (new_kept_tokens, dropped_reason_map) with collapse reasons added.
    """
    by_base: dict[str, list[tuple[int, str, _TailParts]]] = {}
    for idx, tok in enumerate(kept_tokens):
        parsed = _parse_trivial_tail(tok)
        if not parsed:
            continue
        base, tail_parts = parsed
        by_base.setdefault(base, []).append((idx, tok, tail_parts))

    if not by_base:
        return kept_tokens, {}

    token_set = set(kept_tokens)
    to_drop: dict[str, str] = {}

    # Decide per base
    for base, members in by_base.items():
        if len(members) == 1 and base not in token_set:
            # Single member family with no base present -> do nothing
            continue

        if base in token_set:
            # Collapse into the base; drop all family members (not the base)
            for _, tok, _ in members:
                to_drop[tok] = f"collapsed_into_base: base={base}"
        else:
            # Choose a single representative by tail order
            members_sorted = sorted(members, key=lambda x: _family_sort_key(x[2], x[0]))
            keep_tok = members_sorted[0][1]
            for _, tok, _ in members_sorted[1:]:
                to_drop[tok] = f"collapsed_family_no_base: base={base}, kept={keep_tok}"

    # Build final kept list in original order
    new_kept = [t for t in kept_tokens if t not in to_drop]
    return new_kept, to_drop


# ---------- internal/public API (list[str] in/out) ----------
def _simplify_gene_list_internal(
    gene_names: list[str],
    *,
    rules: Rules | None = None,
    with_audit: bool = False,
    collapse_families: bool = True,
) -> tuple[list[str], dict[str, str]]:
    """Internal implementation for gene list simplification.

    Args:
        gene_names: List of raw gene symbols/IDs.
        rules: Optional Rules (conservative by default).
        with_audit: If True, include dropped reasons in return value.
        collapse_families: If True, apply the family collapsing pass.

    Returns:
        (kept_list, dropped_reason_map). When with_audit=False, dropped_reason_map is {}.
    """
    rules = rules or Rules()
    tokens: list[str] = list(
        OrderedDict((t.strip(), None) for t in gene_names if t.strip()).keys()
    )

    # Stage 1: symbol cleanup
    decisions = _token_decisions(tokens, rules)
    if not any(d.kept for d in decisions):  # all uninformative -> keep all
        kept_tokens = [d.token for d in decisions]
        dropped: dict[str, str] = {}
    else:
        kept_tokens = [d.token for d in decisions if d.kept]
        dropped = (
            {d.token: d.reason for d in decisions if not d.kept} if with_audit else {}
        )

    # Stage 2: family collapsing
    if collapse_families and kept_tokens:
        kept_tokens, fam_dropped = _collapse_families(kept_tokens)
        if with_audit and fam_dropped:
            dropped.update(fam_dropped)

    return kept_tokens, dropped


def simplify_gene_list(
    gene_names: list[str],
    *,
    rules: Rules | None = None,
    collapse_families: bool = True,
) -> list[str]:
    """Simplify a list of gene symbols/IDs by removing uninformative tokens and collapsing families.

    Args:
        gene_names: List of gene symbols/IDs.
        rules: Optional Rules (conservative by default).
        collapse_families: If True, apply the family collapsing pass.

    Returns:
        List of kept tokens in first-seen order (post-collapsing).
    """
    kept, _ = _simplify_gene_list_internal(
        gene_names, rules=rules, with_audit=False, collapse_families=collapse_families
    )
    return kept


def simplify_gene_list_with_audit(
    gene_names: list[str],
    *,
    rules: Rules | None = None,
    collapse_families: bool = True,
) -> tuple[list[str], dict[str, str]]:
    """Simplify a list of gene symbols/IDs and also return a map of dropped tokens → reason.

    Args:
        gene_names: List of gene symbols/IDs.
        rules: Optional Rules (conservative by default).
        collapse_families: If True, apply the family collapsing pass.

    Returns:
        (kept_list, dropped_reason_map).
    """
    return _simplify_gene_list_internal(
        gene_names, rules=rules, with_audit=True, collapse_families=collapse_families
    )


# ---------- tests / examples ----------
if __name__ == "__main__":
    conservative = Rules(aggressive_human_p_suffix=False)
    aggressive = Rules(aggressive_human_p_suffix=True)

    def expect(
        inp: list[str],
        out: list[str],
        *,
        rules: Rules = conservative,
        collapse_families: bool = True,
    ) -> None:
        got = simplify_gene_list(inp, rules=rules, collapse_families=collapse_families)
        assert got == out, f"{inp} → {got}, expected {out}"

    # Stage-1 sanity
    expect(["KLRC4-KLRK1", "KLRK1"], ["KLRK1"])
    expect(["IFNAR2-IL10RB", "IL10RB"], ["IL10RB"])
    expect(["IFNAR2-IL10RB"], ["IFNAR2-IL10RB"])
    # all uninformative -> keep all
    expect(["ENSG00000173366", "Gm3839"], ["ENSG00000173366", "Gm3839"])

    # lnc-like / models
    expect(["Gm10358", "Gapdh", "Gm3839"], ["Gapdh"])
    expect(["MMP24", "MMP24OS"], ["MMP24"])
    expect(["FOXP3-AS1", "FOXP3"], ["FOXP3"])

    # Mitochondrial
    expect(["mt-Co3", "Gm3839"], ["mt-Co3"])
    # Look like a gene id but too short, so not treated as gene id and kept
    expect(["MT-CO3", "ENSG0000"], ["MT-CO3", "ENSG0000"])

    # Human-safe pseudogenes (digits before 'P') are dropped.
    expect(["DEFA10P", "DEFA1"], ["DEFA1"])
    expect(["CLCA3P", "CLCA3"], ["CLCA3"])
    expect(["RNA18SP5", "RNA18S"], ["RNA18S"])

    # Aggressive human '...P<digits>' behavior (use carefully!)
    expect(["HMGB1P1", "HMGB1"], ["HMGB1P1", "HMGB1"], rules=conservative)
    expect(["HMGB1P1", "HMGB1"], ["HMGB1"], rules=aggressive)
    expect(["FOXP1", "FOXP3"], ["FOXP1"], rules=conservative)  # not dropped alone
    expect(["FOXP1", "FOXP3"], ["FOXP1"], rules=aggressive)  # not dropped alone
    # Both are pseudogenes, but have a common prefix, so keep the first one
    expect(["KRT8P1", "KRT8P2"], ["KRT8P1"])
    expect(["KRT8P1", "KRT8"], ["KRT8"], rules=aggressive)

    # Histone paralogs remain (no trivial tail)
    expect(["H2BC12", "H2BC12L"], ["H2BC12", "H2BC12L"])

    expect(["Ppp1r1b", "1700003D09Rik"], ["Ppp1r1b", "1700003D09Rik"])

    # ---------- NEW: family collapsing ----------
    # Case A: base present -> collapse to base
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

    # Case B: base absent -> keep earliest by tail order
    expect(["Selenbp1", "Selenbp2"], ["Selenbp1"])

    # Mixed families + outsiders
    # Klk1-family picks 1st; Klk2 untouched
    expect(["Klk1b1", "Klk1b2", "Klk2"], ["Klk1b1", "Klk2"])

    # Do not over-collapse when only 1 trivial-tail member exists
    expect(["Klk1b1", "Gm3839"], ["Klk1b1"])  # single-member family -> keep as-is

    # Do not mix with nontrivial tails (uppercase letters, here an L)
    expect(["H2BC12", "H2BC12L", "H2BC13"], ["H2BC12", "H2BC12L"])

    # Audit demo
    demo = [
        "Klk1",
        "Klk1b1",
        "Klk1b3",
        "Selenbp1",
        "Selenbp2",
        "HMGB1P1",
        "HMGB1",
        "MMP24OS",
        "MMP24",
    ]
    kept, dropped = simplify_gene_list_with_audit(demo, rules=conservative)
    assert kept == ["Klk1", "Selenbp1", "HMGB1P1", "HMGB1", "MMP24"], (kept, dropped)

    print("All tests passed ✓")
