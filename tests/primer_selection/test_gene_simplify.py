from folitools.primer_selection.recover_utils import simplify_gene_list, simplify_gene_list_with_audit


def test_simplify_gene_list_basic_cases():
    """Test basic gene list simplification cases."""
    cases = [
        (["KLRC4-KLRK1", "KLRK1"], ["KLRK1"]),
        (["IFNAR2-IL10RB", "IL10RB"], ["IL10RB"]),
        (["IFNAR2-IL10RB"], ["IFNAR2-IL10RB"]),
        (["ENSG00000173366", "Gm3839"], ["ENSG00000173366", "Gm3839"]),
        (["Gm10358", "Gapdh", "Gm3839"], ["Gapdh"]),
        (["HMGB1", "HMGB1P1"], ["HMGB1", "HMGB1P1"]),  # Conservative mode keeps both
        (["HMGB1P1"], ["HMGB1P1"]),
        (["MMP24", "MMP24OS"], ["MMP24"]),
        (["MMP24OS"], ["MMP24OS"]),
        (["H2BC12", "H2BC12L"], ["H2BC12", "H2BC12L"]),
        # DEFA3 gets collapsed because DEFA + "3" is seen as family member, keeping DEFA1 (earliest)
        (["DEFA1", "DEFA3", "DEFA1B"], ["DEFA1", "DEFA1B"]),
        (["FOXP3-AS1", "FOXP3"], ["FOXP3"]),
        (["KLRC4", "KLRC4-KLRK1", "KLRK1"], ["KLRC4", "KLRK1"]),
        (["KLRK1", "KLRK1"], ["KLRK1"]),
        (["DEFA10P", "DEFA1"], ["DEFA1"]),  # DEFA10P has digits before P, so it's dropped
    ]

    for input_list, expected in cases:
        result = simplify_gene_list(input_list)
        assert result == expected, f"Input: {input_list} → Got: {result!r}, Expected: {expected!r}"


def test_simplify_gene_list_with_audit():
    """Test gene list simplification with audit functionality."""
    demo = ["ENSG00000173366", "Gm3839", "IFNAR2-IL10RB", "IL10RB", "HMGB1P1", "HMGB1", "MMP24OS", "MMP24"]
    simplified, dropped = simplify_gene_list_with_audit(demo)
    
    # Check that we get a simplified result and a dropped dictionary
    assert isinstance(simplified, list)
    assert isinstance(dropped, dict)
    
    # The simplified result should contain informative genes
    assert "HMGB1" in simplified
    assert "MMP24" in simplified
    
    # Check that some tokens were dropped with reasons
    assert len(dropped) > 0
    for token, reason in dropped.items():
        assert isinstance(reason, str)
        assert reason.startswith("dropped:")


def test_empty_and_edge_cases():
    """Test edge cases like empty strings and whitespace."""
    # Empty list
    assert simplify_gene_list([]) == []
    
    # Empty strings in list
    assert simplify_gene_list(["", "  ", ""]) == []
    
    # Whitespace handling - family collapsing applies to GENE1/GENE2
    assert simplify_gene_list([" GENE1 ", " GENE2 "]) == ["GENE1"]  # Family collapsing: GENE base, keep GENE1
    
    # Single gene
    assert simplify_gene_list(["GENE1"]) == ["GENE1"]


def test_readthrough_gene_logic():
    """Test specific readthrough gene behavior."""
    # Readthrough should be dropped when component is present
    assert simplify_gene_list(["KLRC4-KLRK1", "KLRK1"]) == ["KLRK1"]
    
    # Readthrough should be kept when no component is present
    assert simplify_gene_list(["IFNAR2-IL10RB"]) == ["IFNAR2-IL10RB"]
    
    # Multiple components present
    assert simplify_gene_list(["KLRC4", "KLRC4-KLRK1", "KLRK1"]) == ["KLRC4", "KLRK1"]


def test_uninformative_gene_categories():
    """Test that uninformative gene categories are properly handled."""
    # All uninformative - should keep everything
    result = simplify_gene_list(["ENSG00000173366", "Gm3839"])
    assert "ENSG00000173366" in result
    assert "Gm3839" in result
    
    # Mix of informative and uninformative - should keep only informative
    result = simplify_gene_list(["ENSG00000173366", "GAPDH", "Gm3839"])
    assert result == ["GAPDH"]
    
    # Pseudogenes with digits before P (safe patterns)
    result = simplify_gene_list(["DEFA1", "DEFA10P"])
    assert result == ["DEFA1"]
    
    # P<digits> patterns (conservative mode keeps both)
    result = simplify_gene_list(["HMGB1", "HMGB1P1"])
    assert result == ["HMGB1", "HMGB1P1"]
    
    # lncRNA-like patterns
    result = simplify_gene_list(["MMP24", "MMP24OS"])
    assert result == ["MMP24"]
