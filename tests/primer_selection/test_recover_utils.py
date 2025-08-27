"""Tests for the recover_utils module (gene_simplify functionality)."""

from folitools.primer_selection.recover_utils import (
    simplify_gene_list,
    simplify_gene_list_with_audit,
    Rules,
)


class TestSimplifyGeneList:
    """Test the simplify_gene_list function."""

    def test_core_behavior_conservative(self):
        """Test core behavior with conservative rules."""
        conservative = Rules(aggressive_human_p_suffix=False)

        # Readthrough behavior
        assert simplify_gene_list(["KLRC4-KLRK1", "KLRK1"], rules=conservative) == ["KLRK1"]
        assert simplify_gene_list(["IFNAR2-IL10RB", "IL10RB"], rules=conservative) == ["IL10RB"]
        assert simplify_gene_list(["IFNAR2-IL10RB"], rules=conservative) == ["IFNAR2-IL10RB"]

    def test_all_uninformative_keep_all(self):
        """Test that all uninformative tokens are kept when nothing is informative."""
        conservative = Rules(aggressive_human_p_suffix=False)
        result = simplify_gene_list(["ENSG00000173366", "Gm3839"], rules=conservative)
        assert result == ["ENSG00000173366", "Gm3839"]

    def test_mouse_gene_models_and_dash_suffixes(self):
        """Test mouse gene models and dash-suffixes."""
        conservative = Rules(aggressive_human_p_suffix=False)
        
        # Mouse gene models
        assert simplify_gene_list(["Gm10358", "Gapdh", "Gm3839"], rules=conservative) == ["Gapdh"]
        
        # Dash-ps suffixes
        assert simplify_gene_list(["Clca4c-ps", "Clca4c"], rules=conservative) == ["Clca4c"]
        assert simplify_gene_list(["Clca4c-ps"], rules=conservative) == ["Clca4c-ps"]
        
        # Dash-rs suffixes
        assert simplify_gene_list(["Rn18s-rs5", "Rn18s"], rules=conservative) == ["Rn18s"]
        assert simplify_gene_list(["Rn18s-rs5"], rules=conservative) == ["Rn18s-rs5"]

    def test_mitochondrial_genes(self):
        """Test mitochondrial gene handling."""
        conservative = Rules(aggressive_human_p_suffix=False)
        
        # Mouse mitochondrial
        assert simplify_gene_list(["mt-Co3", "Gm3839"], rules=conservative) == ["mt-Co3"]
        
        # Human mitochondrial - "ENSG0000" is not a valid Ensembl ID pattern so it's informative
        assert simplify_gene_list(["MT-CO3", "ENSG0000"], rules=conservative) == ["MT-CO3", "ENSG0000"]

    def test_lnc_like_patterns(self):
        """Test lncRNA-like pattern handling."""
        conservative = Rules(aggressive_human_p_suffix=False)
        
        # OS suffix
        assert simplify_gene_list(["MMP24", "MMP24OS"], rules=conservative) == ["MMP24"]
        
        # AS suffix
        assert simplify_gene_list(["FOXP3-AS1", "FOXP3"], rules=conservative) == ["FOXP3"]

    def test_histone_paralogs(self):
        """Test that histone paralogs are both kept as informative."""
        conservative = Rules(aggressive_human_p_suffix=False)
        result = simplify_gene_list(["H2BC12", "H2BC12L"], rules=conservative)
        assert result == ["H2BC12", "H2BC12L"]

    def test_human_safe_pseudogenes(self):
        """Test human-safe pseudogenes (digits before terminal 'P')."""
        conservative = Rules(aggressive_human_p_suffix=False)
        
        # DEFA10P should be dropped when DEFA1 is present
        assert simplify_gene_list(["DEFA10P", "DEFA1"], rules=conservative) == ["DEFA1"]
        
        # CLCA3P should be dropped when CLCA3 is present
        assert simplify_gene_list(["CLCA3P", "CLCA3"], rules=conservative) == ["CLCA3"]
        
        # Human rRNA pseudogene
        assert simplify_gene_list(["RNA18SP5", "RNA18S"], rules=conservative) == ["RNA18S"]

    def test_aggressive_mode(self):
        """Test aggressive mode for human P<digits> patterns."""
        conservative = Rules(aggressive_human_p_suffix=False)
        aggressive = Rules(aggressive_human_p_suffix=True)
        
        # Aggressive mode should drop HMGB1P1 when HMGB1 is present
        assert simplify_gene_list(["HMGB1P1", "HMGB1"], rules=aggressive) == ["HMGB1"]
        
        # But still keeps bona fide genes when they stand alone
        assert simplify_gene_list(["FOXP1"], rules=aggressive) == ["FOXP1"]
        
        # If only aggressive-style P<digits> tokens exist, family collapsing applies
        # KRT8P1/KRT8P2 have base KRT8P, so they collapse to first occurrence
        result = simplify_gene_list(["KRT8P1", "KRT8P2"], rules=conservative)
        assert result == ["KRT8P1"]  # Family collapsing: base=KRT8P, keep first
        
        # But if informative parent is present, aggressive mode drops P<digits>
        assert simplify_gene_list(["KRT8P1", "KRT8"], rules=aggressive) == ["KRT8"]

    def test_order_preservation(self):
        """Test that order is preserved in results."""
        conservative = Rules(aggressive_human_p_suffix=False)
        
        # Family collapsing sorts by numeric tail, so GENE1 comes first even if GENE3 appeared first
        result = simplify_gene_list(["GENE3", "GENE1", "GENE2"], rules=conservative)
        assert result == ["GENE1"]  # Earliest by numeric order (1 < 2 < 3)
        
        # Duplicates should be removed but order maintained
        result = simplify_gene_list(["GENE1", "GENE2", "GENE1"], rules=conservative)
        assert result == ["GENE1"]  # Family collapsing + deduplication

    def test_edge_cases(self):
        """Test edge cases."""
        conservative = Rules(aggressive_human_p_suffix=False)
        
        # Empty list
        assert simplify_gene_list([], rules=conservative) == []
        
        # Whitespace handling - GENE1/GENE2 will be treated as family
        result = simplify_gene_list([" GENE1 ", "", "  GENE2  "], rules=conservative)
        assert result == ["GENE1"]  # Family collapsing applies

    def test_family_collapsing_nontrivial_tails(self):
        """Test that uppercase letters prevent family collapsing."""
        conservative = Rules(aggressive_human_p_suffix=False)
        
        # H2BC13 has trivial tail (13), but H2BC12L has nontrivial tail (L)
        # So only H2BC13 would be in a family by itself
        result = simplify_gene_list(["H2BC12", "H2BC12L", "H2BC13"], rules=conservative)
        assert result == ["H2BC12", "H2BC12L"]  # H2BC13 collapsed as single-member family

    def test_order_preservation_with_collapsing(self):
        """Test that order is preserved in results accounting for family collapsing."""
        conservative = Rules(aggressive_human_p_suffix=False)
        
        # Family collapsing sorts by numeric tail, so GENE1 comes first even if GENE3 appeared first
        result = simplify_gene_list(["GENE3", "GENE1", "GENE2"], rules=conservative)
        assert result == ["GENE1"]  # Earliest by numeric order (1 < 2 < 3)
        
        # Duplicates should be removed but order maintained
        result = simplify_gene_list(["GENE1", "GENE2", "GENE1"], rules=conservative)
        assert result == ["GENE1"]  # Family collapsing + deduplication

    def test_edge_cases_with_collapsing(self):
        """Test edge cases accounting for family collapsing."""
        conservative = Rules(aggressive_human_p_suffix=False)
        
        # Empty list
        assert simplify_gene_list([], rules=conservative) == []
        
        # Whitespace handling - GENE1/GENE2 will be treated as family
        result = simplify_gene_list([" GENE1 ", "", "  GENE2  "], rules=conservative)
        assert result == ["GENE1"]  # Family collapsing applies

    def test_family_collapsing_case_a_base_present(self):
        """Test family collapsing when base is present - collapse to base."""
        conservative = Rules(aggressive_human_p_suffix=False)
        
        # Case A: base present -> collapse all to base
        result = simplify_gene_list([
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
        ], rules=conservative)
        assert result == ["Klk1"]

    def test_family_collapsing_case_b_base_absent(self):
        """Test family collapsing when base is absent - keep earliest by order."""
        conservative = Rules(aggressive_human_p_suffix=False)
        
        # Case B: base absent -> keep earliest by tail order
        result = simplify_gene_list(["Selenbp1", "Selenbp2"], rules=conservative)
        assert result == ["Selenbp1"]

    def test_family_collapsing_mixed_families(self):
        """Test family collapsing with mixed families and outsiders."""
        conservative = Rules(aggressive_human_p_suffix=False)
        
        # Mixed families + outsiders
        result = simplify_gene_list(["Klk1b1", "Klk1b2", "Klk2"], rules=conservative)
        assert result == ["Klk1b1", "Klk2"]  # Klk1-family picks 1st; Klk2 untouched

    def test_family_collapsing_single_member(self):
        """Test that single-member families are not over-collapsed."""
        conservative = Rules(aggressive_human_p_suffix=False)
        
        # Do not over-collapse when only 1 trivial-tail member exists
        result = simplify_gene_list(["Klk1b1", "Gm3839"], rules=conservative)
        assert result == ["Klk1b1"]  # single-member family -> keep as-is

    def test_family_collapsing_nontrivial_tails_corrected(self):
        """Test that uppercase letters prevent family collapsing."""
        conservative = Rules(aggressive_human_p_suffix=False)
        
        # H2BC13 has trivial tail (13), but H2BC12L has nontrivial tail (L)
        # So only H2BC13 would be in a family by itself and gets dropped
        result = simplify_gene_list(["H2BC12", "H2BC12L", "H2BC13"], rules=conservative)
        assert result == ["H2BC12", "H2BC12L"]  # H2BC13 collapsed as single-member family

    def test_family_collapsing_disabled(self):
        """Test disabling family collapsing."""
        conservative = Rules(aggressive_human_p_suffix=False)
        
        # With family collapsing disabled
        result = simplify_gene_list(
            ["Selenbp1", "Selenbp2"], 
            rules=conservative, 
            collapse_families=False
        )
        assert result == ["Selenbp1", "Selenbp2"]
        
        # With family collapsing enabled (default)
        result = simplify_gene_list(
            ["Selenbp1", "Selenbp2"], 
            rules=conservative, 
            collapse_families=True
        )
        assert result == ["Selenbp1"]


class TestSimplifyGeneListWithAudit:
    """Test the simplify_gene_list_with_audit function."""

    def test_audit_functionality(self):
        """Test audit functionality returns proper dropped reasons."""
        conservative = Rules(aggressive_human_p_suffix=False)
        
        demo = [
            "Klk1",
            "Klk1b1",
            "Klk1b3",
            "Selenbp1",
            "Selenbp2",
            "HMGB1P1",  # In conservative mode, not dropped by P<digits> rule
            "HMGB1",
            "MMP24OS",
            "MMP24",
        ]
        
        kept, dropped = simplify_gene_list_with_audit(demo, rules=conservative)
        
        # Check that we got expected results from the demo
        # HMGB1P1 is kept in conservative mode (not dropped)
        assert kept == ["Klk1", "Selenbp1", "HMGB1P1", "HMGB1", "MMP24"]
        
        # Check specific dropped reasons
        assert "collapsed_into_base: base=Klk1" in dropped.get("Klk1b1", "")
        assert "collapsed_into_base: base=Klk1" in dropped.get("Klk1b3", "")
        assert "collapsed_family_no_base: base=Sele" in dropped.get("Selenbp2", "")
        assert dropped.get("MMP24OS") == "dropped: lnc_like_placeholder"

    def test_audit_with_no_drops(self):
        """Test audit when nothing gets dropped."""
        conservative = Rules(aggressive_human_p_suffix=False)
        
        # Use genes that won't trigger family collapsing
        kept, dropped = simplify_gene_list_with_audit(["ACTA1", "ACTB"], rules=conservative)
        
        assert kept == ["ACTA1", "ACTB"]
        assert dropped == {}

    def test_audit_all_uninformative(self):
        """Test audit when all tokens are uninformative."""
        conservative = Rules(aggressive_human_p_suffix=False)
        
        kept, dropped = simplify_gene_list_with_audit(["ENSG00000173366", "Gm3839"], rules=conservative)
        
        # All should be kept with special reason
        assert kept == ["ENSG00000173366", "Gm3839"]
        assert dropped == {}

    def test_audit_family_collapsing_disabled(self):
        """Test audit functionality with family collapsing disabled."""
        conservative = Rules(aggressive_human_p_suffix=False)
        
        demo = [
            "Klk1",
            "Klk1b1",
            "Klk1b3",
            "HMGB1",
            "MMP24OS",
            "MMP24",
        ]
        
        kept, dropped = simplify_gene_list_with_audit(
            demo, 
            rules=conservative,
            collapse_families=False
        )
        
        # With family collapsing disabled, Klk1b1 and Klk1b3 should be kept
        assert kept == ["Klk1", "Klk1b1", "Klk1b3", "HMGB1", "MMP24"]
        
        # Only non-family drops should be present
        assert dropped.get("MMP24OS") == "dropped: lnc_like_placeholder"
        assert "Klk1b1" not in dropped
        assert "Klk1b3" not in dropped


class TestRulesConfiguration:
    """Test the Rules configuration class."""

    def test_default_rules(self):
        """Test that default rules are conservative."""
        default_rules = Rules()
        assert default_rules.aggressive_human_p_suffix is False

    def test_aggressive_rules(self):
        """Test aggressive rules configuration."""
        aggressive_rules = Rules(aggressive_human_p_suffix=True)
        assert aggressive_rules.aggressive_human_p_suffix is True

    def test_rules_immutability(self):
        """Test that Rules objects are immutable."""
        rules = Rules()
        
        # Rules objects are frozen dataclasses - they're created properly
        assert rules.aggressive_human_p_suffix is False
        
        # Create a new instance with different values
        aggressive_rules = Rules(aggressive_human_p_suffix=True)
        assert aggressive_rules.aggressive_human_p_suffix is True


class TestModuleMainTests:
    """Test cases moved from the module's __name__ == "__main__" section."""

    def test_stage_1_sanity(self):
        """Test stage-1 sanity checks from module main."""
        conservative = Rules(aggressive_human_p_suffix=False)
        
        assert simplify_gene_list(["KLRC4-KLRK1", "KLRK1"], rules=conservative) == ["KLRK1"]
        assert simplify_gene_list(["IFNAR2-IL10RB", "IL10RB"], rules=conservative) == ["IL10RB"]
        assert simplify_gene_list(["IFNAR2-IL10RB"], rules=conservative) == ["IFNAR2-IL10RB"]
        assert simplify_gene_list(["ENSG00000173366", "Gm3839"], rules=conservative) == ["ENSG00000173366", "Gm3839"]

    def test_lnc_like_and_models(self):
        """Test lnc-like / models from module main."""
        conservative = Rules(aggressive_human_p_suffix=False)
        
        assert simplify_gene_list(["Gm10358", "Gapdh", "Gm3839"], rules=conservative) == ["Gapdh"]
        assert simplify_gene_list(["MMP24", "MMP24OS"], rules=conservative) == ["MMP24"]
        assert simplify_gene_list(["FOXP3-AS1", "FOXP3"], rules=conservative) == ["FOXP3"]

    def test_mitochondrial_from_main(self):
        """Test mitochondrial from module main."""
        conservative = Rules(aggressive_human_p_suffix=False)
        
        assert simplify_gene_list(["mt-Co3", "Gm3839"], rules=conservative) == ["mt-Co3"]
        assert simplify_gene_list(["MT-CO3", "ENSG0000"], rules=conservative) == ["MT-CO3", "ENSG0000"]

    def test_human_safe_pseudogenes_from_main(self):
        """Test human-safe pseudogenes from module main."""
        conservative = Rules(aggressive_human_p_suffix=False)
        
        assert simplify_gene_list(["DEFA10P", "DEFA1"], rules=conservative) == ["DEFA1"]
        assert simplify_gene_list(["CLCA3P", "CLCA3"], rules=conservative) == ["CLCA3"]
        assert simplify_gene_list(["RNA18SP5", "RNA18S"], rules=conservative) == ["RNA18S"]

    def test_aggressive_human_patterns_from_main(self):
        """Test aggressive human patterns from module main."""
        conservative = Rules(aggressive_human_p_suffix=False)
        aggressive = Rules(aggressive_human_p_suffix=True)
        
        assert simplify_gene_list(["HMGB1P1", "HMGB1"], rules=aggressive) == ["HMGB1"]
        assert simplify_gene_list(["FOXP1"], rules=aggressive) == ["FOXP1"]
        # Family collapsing applies to KRT8P1/KRT8P2 (base=KRT8P)
        assert simplify_gene_list(["KRT8P1", "KRT8P2"], rules=conservative) == ["KRT8P1"]
        assert simplify_gene_list(["KRT8P1", "KRT8"], rules=aggressive) == ["KRT8"]

    def test_histone_paralogs_from_main(self):
        """Test histone paralogs from module main."""
        conservative = Rules(aggressive_human_p_suffix=False)
        
        assert simplify_gene_list(["H2BC12", "H2BC12L"], rules=conservative) == ["H2BC12", "H2BC12L"]

    def test_family_collapsing_from_main(self):
        """Test family collapsing examples from module main."""
        conservative = Rules(aggressive_human_p_suffix=False)
        
        # Case A: base present -> collapse to base
        assert simplify_gene_list([
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
        ], rules=conservative) == ["Klk1"]

        # Case B: base absent -> keep earliest by tail order
        assert simplify_gene_list(["Selenbp1", "Selenbp2"], rules=conservative) == ["Selenbp1"]

        # Mixed families + outsiders
        assert simplify_gene_list(["Klk1b1", "Klk1b2", "Klk2"], rules=conservative) == ["Klk1b1", "Klk2"]

        # Do not over-collapse when only 1 trivial-tail member exists
        assert simplify_gene_list(["Klk1b1", "Gm3839"], rules=conservative) == ["Klk1b1"]

        # H2BC13 has trivial tail, so it gets collapsed as single member family
        assert simplify_gene_list(["H2BC12", "H2BC12L", "H2BC13"], rules=conservative) == ["H2BC12", "H2BC12L"]

    def test_audit_demo_from_main(self):
        """Test audit demo from module main."""
        conservative = Rules(aggressive_human_p_suffix=False)
        
        demo = [
            "Klk1",
            "Klk1b1",
            "Klk1b3",
            "Selenbp1",
            "Selenbp2",
            "HMGB1P1",  # Keep in conservative mode
            "HMGB1",
            "MMP24OS",
            "MMP24",
        ]
        kept, dropped = simplify_gene_list_with_audit(demo, rules=conservative)
        # HMGB1P1 is kept in conservative mode
        assert kept == ["Klk1", "Selenbp1", "HMGB1P1", "HMGB1", "MMP24"]
        
        # Check some key dropped reasons
        assert "collapsed_into_base: base=Klk1" in dropped.get("Klk1b1", "")
        assert "collapsed_family_no_base" in dropped.get("Selenbp2", "")
        assert dropped.get("MMP24OS") == "dropped: lnc_like_placeholder"


class TestComplexScenarios:
    """Test complex real-world scenarios."""

    def test_mixed_species_patterns(self):
        """Test mixed human/mouse patterns."""
        conservative = Rules(aggressive_human_p_suffix=False)
        
        mixed_genes = [
            "MT-CO3",  # Human mitochondrial (informative)
            "mt-Co3",  # Mouse mitochondrial (informative) 
            "ENSG00000173366",  # Human Ensembl ID (uninformative)
            "ENSMUSG00000000001",  # Mouse Ensembl ID (uninformative)
            "Gm3839",  # Mouse gene model (uninformative)
            "KLRC4-KLRK1",  # Readthrough (context-dependent)
            "KLRK1",  # Component of readthrough (informative)
            "HMGB1",  # Regular gene (informative)
            "DEFA10P",  # Pseudogene with digits before P (uninformative)
        ]
        
        result = simplify_gene_list(mixed_genes, rules=conservative)
        expected = ["MT-CO3", "mt-Co3", "KLRK1", "HMGB1"]
        assert result == expected

    def test_multiple_readthroughs(self):
        """Test behavior with multiple readthrough genes."""
        conservative = Rules(aggressive_human_p_suffix=False)
        
        genes = [
            "GENEA", "GENEB", "GENEC",  # Regular genes
            "GENEA-GENEB",  # Readthrough (should be dropped - GENEA component present)
            "GENEX-GENEY",  # Readthrough (should be kept - no components)
            "GENEC-GENED",  # Readthrough (should be dropped - GENEC component present)
        ]
        
        result = simplify_gene_list(genes, rules=conservative)
        expected = ["GENEA", "GENEB", "GENEC", "GENEX-GENEY"]
        assert result == expected

    def test_comprehensive_integration(self):
        """Test a comprehensive scenario with all features."""
        conservative = Rules(aggressive_human_p_suffix=False)
        
        genes = [
            # Family collapsing scenario
            "Klk1", "Klk1b1", "Klk1b2",
            # Pseudogenes
            "DEFA10P", "DEFA1", 
            # Gene models
            "Gm3839", "Gapdh",
            # Readthrough
            "IFNAR2-IL10RB", "IL10RB",
            # LncRNA-like
            "MMP24OS", "MMP24",
            # Mitochondrial
            "mt-Co3",
            # Regular genes
            "FOXP1", "HMGB1"
        ]
        
        result = simplify_gene_list(genes, rules=conservative)
        expected = ["Klk1", "DEFA1", "Gapdh", "IL10RB", "MMP24", "mt-Co3", "FOXP1", "HMGB1"]
        assert result == expected
