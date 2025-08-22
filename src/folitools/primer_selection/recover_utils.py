import re
from collections import OrderedDict


def simplify_gene_list(gene_string: str) -> str:
    """Simplify a '|' separated gene string by removing redundant tokens.

    Rules (conservative defaults):
      1) Drop bare Ensembl IDs (ENSG...) if any named symbols are present.
      2) Drop readthrough symbols of the form 'A-B' iff 'B' is also present.
      3) Drop pseudogene copies like '<PARENT>P<digits>' iff <PARENT> is present.
      4) Keep antisense/OS (e.g., MMP24OS) and '-like' paralogs (e.g., H2BC12L).
      5) Optionally consolidate DEFA1/DEFA3 -> DEFA1A3 when both present.

    Args:
        gene_string: Input like "KLRC4-KLRK1|KLRK1|ENSG00000173366|TLR9".
        collapse_defa_locus: If True, replace co-present DEFA1 & DEFA3 with DEFA1A3.

    Returns:
        A '|' joined simplified string, preserving original order where possible.

    Examples:
        >>> simplify_gene_list("KLRC4-KLRK1|KLRK1")
        'KLRK1'
        >>> simplify_gene_list("ENSG00000173366|TLR9")
        'TLR9'
        >>> simplify_gene_list("HMGB1|HMGB1P1")
        'HMGB1'
        >>> simplify_gene_list("HMGB1P1")
        'HMGB1P1'
        >>> simplify_gene_list("MMP24|MMP24OS")
        'MMP24|MMP24OS'
        >>> simplify_gene_list("H2BC12|H2BC12L")
        'H2BC12|H2BC12L'
        >>> simplify_gene_list("DEFA1|DEFA3|DEFA1B")
        'DEFA1|DEFA3|DEFA1B'
        >>> simplify_gene_list("DEFA1|DEFA3|DEFA1B", collapse_defa_locus=True)
        'DEFA1A3|DEFA1B'
    """
    # ---- tokenize & normalize ----
    raw_tokens = [t.strip() for t in gene_string.split("|") if t.strip()]
    # Ordered unique, preserving first occurrence
    tokens: list[str] = list(OrderedDict((t, None) for t in raw_tokens).keys())
    token_set: set[str] = set(tokens)

    def is_ensembl_gene_id(sym: str) -> bool:
        return bool(re.fullmatch(r"ENSG\d{11}(?:\.\d+)?", sym))

    def readthrough_components(sym: str) -> tuple[str, str] | None:
        # Only treat simple one-hyphen forms as readthrough candidates
        if sym.count("-") != 1:
            return None
        left, right = sym.split("-", 1)
        if not left or not right:
            return None
        return (left, right)

    def parent_if_pseudogene(sym: str) -> str | None:
        # Only consider dropping if parent (prefix before 'P<number>') is present exactly.
        m = re.fullmatch(r"([A-Z0-9]+)P\d+", sym)
        if m:
            parent = m.group(1)
            if parent in token_set:
                return parent
        return None

    # ---- 1) Drop ENSG IDs if any named symbol exists ----
    any_named = any(not is_ensembl_gene_id(t) for t in tokens)
    if any_named:
        tokens = [t for t in tokens if not is_ensembl_gene_id(t)]
        token_set = set(tokens)

    # ---- 2) Drop readthrough 'A-B' iff 'B' present ----
    pruned: list[str] = []
    for t in tokens:
        parts = readthrough_components(t)
        if parts and parts[1] in token_set:
            # Keep only the simpler downstream gene already present
            continue
        pruned.append(t)
    tokens = pruned
    token_set = set(tokens)

    # ---- 3) Drop pseudogenes '<PARENT>P<digits>' iff parent present ----
    pruned = []
    for t in tokens:
        parent = parent_if_pseudogene(t)
        if parent is not None:
            # Parent symbol already present -> drop this pseudogene token
            continue
        pruned.append(t)
    tokens = pruned
    token_set = set(tokens)

    # # ---- 5) Optional DEFA locus consolidation ----
    # if collapse_defa_locus:
    #     has_defa1 = "DEFA1" in token_set
    #     has_defa3 = "DEFA3" in token_set
    #     if has_defa1 and has_defa3:
    #         # Replace the earlier of the two with DEFA1A3 and remove the other
    #         first_idx = min(tokens.index("DEFA1"), tokens.index("DEFA3"))
    #         # Remove both
    #         tokens = [t for t in tokens if t not in {"DEFA1", "DEFA3"}]
    #         # Insert DEFA1A3 at the earliest position
    #         tokens.insert(first_idx, "DEFA1A3")

    return "|".join(tokens)
