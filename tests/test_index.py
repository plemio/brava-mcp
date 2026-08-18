"""Entity resolution and the bundled metadata."""

import pytest

from brava import index as ix
from brava.constants import ANCESTRIES, MAFS, MASKS


class TestGeneResolution:
    @pytest.mark.parametrize("q", ["PCSK9", "pcsk9", "  PCSK9  ", "ENSG00000169174", "ensg00000169174"])
    def test_accepts_what_a_human_types(self, q):
        assert ix.gene_info(ix.resolve_gene(q))["ensg"] == "ENSG00000169174"

    def test_unknown_gene_is_none_not_a_wrong_answer(self):
        assert ix.resolve_gene("NOTAGENE") is None
        assert ix.resolve_gene("") is None

    def test_exact_symbol_outranks_a_substring_coincidence(self):
        # Without explicit ranking a substring hit can win purely by sitting
        # earlier in the array — upstream hit and fixed exactly this.
        hits = ix.search_genes("PCSK9", 10)
        assert hits[0]["gene"] == "PCSK9"
        assert hits[0]["match"] == "exact symbol"

    def test_prefix_search_returns_the_family(self):
        symbols = {h["gene"] for h in ix.search_genes("PCSK", 10)}
        assert {"PCSK9", "PCSK1"} <= symbols

    def test_gene_info_carries_grch38_coordinates(self):
        info = ix.gene_info(ix.resolve_gene("PCSK9"))
        assert info["chr"] == "1"
        assert info["start"] < info["end"]


class TestPhenotypeResolution:
    @pytest.mark.parametrize("q", ["LDLC", "ldlc", "LDL cholesterol"])
    def test_accepts_id_or_name(self, q):
        assert ix.phenotypes()[ix.resolve_phenotype(q)]["id"] == "LDLC"

    def test_ambiguous_substring_refuses_to_guess(self):
        # "cholesterol" matches LDL, HDL and total — answering one of them
        # silently would be worse than saying nothing.
        assert ix.resolve_phenotype("cholesterol") is None

    def test_search_surfaces_the_candidates_instead(self):
        hits = ix.search_phenotypes("cholesterol", 10)
        assert {h["trait_id"] for h in hits} >= {"LDLC", "HDLC", "TChol"}

    def test_catalogue_has_the_44_published_traits(self):
        assert len(ix.phenotypes()) == 44


class TestVocabularyResolution:
    @pytest.mark.parametrize("value,expected", [("pLoF", 0), ("all variants", 5), (MASKS[4], 4)])
    def test_masks_accept_label_or_raw_string(self, value, expected):
        assert ix.resolve_mask(value) == expected

    @pytest.mark.parametrize("value,expected", [(0.001, 0), ("<0.01%", 1), ("0.0001", 1), (MAFS[0], 0)])
    def test_maf_accepts_number_or_label(self, value, expected):
        assert ix.resolve_maf(value) == expected

    @pytest.mark.parametrize("value,expected", [("skato", "SKAT-O"), ("SKAT-O", "SKAT-O"), ("burden", "Burden")])
    def test_tests_are_case_and_hyphen_insensitive(self, value, expected):
        assert ix.resolve_test(value, "SKAT-O") == expected

    def test_none_means_no_filter_not_an_error(self):
        assert ix.resolve_mask(None) is None
        assert ix.resolve_maf(None) is None
        assert ix.resolve_ancestry(None) is None

    @pytest.mark.parametrize("bad", ["nonsense", "pLOFF"])
    def test_unknown_values_name_the_valid_ones(self, bad):
        with pytest.raises(ValueError) as exc:
            ix.resolve_mask(bad)
        assert "Valid" in str(exc.value)

    def test_non_eur_spelling_variants(self):
        assert ix.resolve_ancestry("non_EUR") == ANCESTRIES.index("non_EUR")
        assert ix.resolve_ancestry("non-EUR") == ANCESTRIES.index("non_EUR")


class TestBundledMetadata:
    def test_every_ancestry_shard_is_present(self):
        for anc in ANCESTRIES:
            assert ix.all_results(anc)["n"] > 0

    def test_biobanks_cover_the_consortium(self):
        assert len(ix.biobanks()) >= 8
        assert all(b["sample_size"] > 0 for b in ix.biobanks())
