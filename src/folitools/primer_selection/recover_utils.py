# gene_simplify.py
from dataclasses import dataclass
from collections import OrderedDict
import re


# ---------- patterns (compiled once) ----------
_ENS_GENE_RE = re.compile(
    r"ENS[A-Z]*G\d{11}(?:\.\d+)?$"
)  # e.g., ENSG..., ENSMUSG..., optional version
_GENE_MODEL_RE = re.compile(
    r"Gm\d+(?:[A-Za-z0-9]+)?$"
)  # mouse "gene model" placeholders: Gm10358, Gm3839
_PSEUDOGENE_SUFFIX_RE = re.compile(
    r"[A-Z0-9]{4,}P\d+$"
)  # HMGB1P1, RPL23AP82, TP73P1 (avoid false hits like MMP24)
_LNC_LIKE_PATTERNS = [
    re.compile(r".*OS\d*$"),  # MMP24OS, ...OS2
    re.compile(r".*-AS\d*$"),  # FOXP3-AS1, ...-AS2
    re.compile(r".*-AS$"),  # ...-AS
    re.compile(r"^LINC\d+[A-Z]*$"),  # LINC01234
    re.compile(
        r"^(AC|AL)\d{3,}(?:\.\d+)?$"
    ),  # AC123456.1 / AL123456.1 (common lnc placeholders)
    re.compile(r"^RP\d+-\d+(?:\.\d+)?$"),  # RP11-... style placeholders
]
_READTHROUGH_SPLIT_RE = re.compile(
    r"([A-Za-z0-9]+)-([A-Za-z0-9]+)$"
)  # simple A-B readthrough symbol


@dataclass(frozen=True)
class TokenDecision:
    """Decision record for a single token."""

    token: str
    kept: bool
    reason: str  # e.g., "informative", "dropped: ensembl_id", "dropped: readthrough_component_present", ...


# ---------- helpers ----------
def _is_ensembl_gene_id(symbol: str) -> bool:
    return bool(_ENS_GENE_RE.fullmatch(symbol))


def _is_gene_model(symbol: str) -> bool:
    return bool(_GENE_MODEL_RE.fullmatch(symbol))


def _is_pseudogene_by_suffix(symbol: str) -> bool:
    return bool(_PSEUDOGENE_SUFFIX_RE.fullmatch(symbol))


def _is_lnc_like_symbol(symbol: str) -> bool:
    return any(p.fullmatch(symbol) for p in _LNC_LIKE_PATTERNS)


def _readthrough_parts(symbol: str) -> tuple[str, str] | None:
    """Return (left, right) if `symbol` looks like a simple readthrough A-B; otherwise None."""
    m = _READTHROUGH_SPLIT_RE.fullmatch(symbol)
    if not m:
        return None
    left, right = m.group(1), m.group(2)
    return (left, right) if left and right else None


def _token_decisions(tokens: list[str]) -> list[TokenDecision]:
    """
    Classify tokens as informative or not, applying context-aware rules for readthroughs.

    Rules (as requested):
      • Pseudogene / lncRNA-like / gene models / plain Ensembl IDs are uninformative in all contexts.
      • Readthrough 'A-B' is uninformative only when any component (A or B) is present;
        otherwise it is considered informative.
      • Everything else is informative.
    """
    token_set = set(tokens)
    decisions: list[TokenDecision] = []

    for tok in tokens:
        # hard uninformative classes
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
        if _is_pseudogene_by_suffix(tok):
            decisions.append(
                TokenDecision(tok, kept=False, reason="dropped: pseudogene_suffix")
            )
            continue
        if _is_lnc_like_symbol(tok):
            decisions.append(
                TokenDecision(tok, kept=False, reason="dropped: lnc_like_placeholder")
            )
            continue

        # readthrough (context-aware)
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

        # default informative
        decisions.append(TokenDecision(tok, kept=True, reason="informative"))

    # if no informative tokens at all, keep everything (mark reason accordingly)
    if not any(d.kept for d in decisions):
        return [
            TokenDecision(d.token, kept=True, reason="kept: all_tokens_uninformative")
            for d in decisions
        ]

    return decisions


# ---------- public API ----------
def _simplify_gene_list_internal(
    gene_string: str, *, with_audit: bool = False
) -> tuple[str, dict[str, str]]:
    """Internal implementation for gene list simplification.

    Args:
        gene_string: Input like "KLRC4-KLRK1|KLRK1|ENSG00000173366|TLR9".
        with_audit: If True, include dropped reasons in return value.

    Returns:
        A tuple (simplified_string, dropped_reason_map).
        When with_audit=False, dropped_reason_map will be empty.
    """
    raw_tokens = [t.strip() for t in gene_string.split("|") if t.strip()]
    tokens: list[str] = list(OrderedDict((t, None) for t in raw_tokens).keys())

    decisions = _token_decisions(tokens)
    if any(d.reason == "kept: all_tokens_uninformative" for d in decisions):
        kept_tokens = [d.token for d in decisions]  # keep all (deduped)
        dropped = {}
    else:
        kept_tokens = [d.token for d in decisions if d.kept]
        dropped = (
            {d.token: d.reason for d in decisions if not d.kept} if with_audit else {}
        )

    simplified_string = "|".join(kept_tokens)
    return simplified_string, dropped


def simplify_gene_list(gene_string: str) -> str:
    """Simplify a '|' separated gene string by removing redundant/uninformative tokens.

    Behavior:
      • If there is at least one informative symbol, keep all and only the informative ones.
      • If all tokens are uninformative, keep everything (order-stable, de-duplicated).
      • Readthrough 'A-B' is uninformative only when any component is present; otherwise informative.
      • Pseudogene / lncRNA-like / gene models / plain Ensembl IDs are uninformative regardless.

    Args:
        gene_string: Input like "KLRC4-KLRK1|KLRK1|ENSG00000173366|TLR9".

    Returns:
        A '|' joined simplified string, preserving first-seen order among kept tokens.
    """
    simplified_string, _ = _simplify_gene_list_internal(gene_string, with_audit=False)
    return simplified_string


def simplify_gene_list_with_audit(gene_string: str) -> tuple[str, dict[str, str]]:
    """Simplify a gene string and also return a map of dropped tokens → reason.

    Args:
        gene_string: Input like "KLRC4-KLRK1|KLRK1|ENSG00000173366|TLR9".

    Returns:
        A pair (simplified_string, dropped_reason_map), where:
          • simplified_string is the '|' joined kept tokens.
          • dropped_reason_map maps each dropped token to its reason.
    """
    return _simplify_gene_list_internal(gene_string, with_audit=True)
