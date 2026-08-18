"""Decoding correctness: the columnar integer encoding is where errors hide."""

import math

import pytest

from brava import index as ix, query as q
from brava.constants import MASK_LABEL, p_from_lp, tier


class TestEffectInterpretation:
    def test_binary_beta_sign_maps_to_risk_or_protection(self):
        assert q.effect_label(0.3, "binary") == "risk-increasing"
        assert q.effect_label(-0.3, "binary") == "protective"

    def test_quantitative_avoids_risk_language(self):
        # A lower LDL-C is not "protective" in the data's own terms: the trait
        # has no polarity, so the label must not invent one.
        assert q.effect_label(-0.037, "quantitative") == "lowers"
        assert q.effect_label(0.037, "quantitative") == "raises"

    def test_zero_and_missing_have_no_direction(self):
        assert q.effect_label(0.0, "binary") is None
        assert q.effect_label(None, "binary") is None

    def test_odds_ratio_only_makes_sense_with_a_ci(self):
        orr = q.odds_ratio(0.3, 0.05)
        assert orr["or"] == pytest.approx(math.exp(0.3), rel=1e-3)
        assert orr["or_lo"] < orr["or"] < orr["or_hi"]
        assert "or_lo" not in q.odds_ratio(0.3, None)


class TestPValueRendering:
    def test_extreme_p_stays_compact(self):
        # A literal decimal here costs ~170 characters per cell.
        assert q.fmt_p(1.584893e-159) == "1.58e-159"
        assert len(q.fmt_p(1.584893e-159)) < 12

    def test_underflow_is_reported_as_a_floor_not_as_zero(self):
        # Upstream floors underflowed p-values rather than nulling them: this is
        # the MOST significant result, and must never read as "no signal".
        assert q.fmt_p(0.0) == "<5e-324"

    def test_missing_stays_missing(self):
        assert q.fmt_p(None) is None

    def test_p_from_lp_round_trips(self):
        assert p_from_lp(158.8) == pytest.approx(1.5849e-159, rel=1e-3)
        assert p_from_lp(None) is None


class TestSignificanceTiers:
    @pytest.mark.parametrize(
        "p,expected",
        [
            (1e-200, "exome-wide"),
            (1.0e-7, "exome-wide"),
            (3e-7, "gene-level"),
            (1e-5, "suggestive"),
            (0.2, "ns"),
            (None, "n/a"),
        ],
    )
    def test_tier_boundaries(self, p, expected):
        assert tier(p) == expected


class TestGeneRows:
    def test_decodes_indices_into_labels(self, pcsk9_gene):
        rows = q.gene_rows(
            pcsk9_gene,
            ix.phenotypes(),
            ancestry_idx=0,
            mask_idx=0,
            maf_idx=0,
            test="SKAT-O",
            max_p=None,
        )
        top = rows[0]
        assert top["trait_id"] == "LDLC"
        assert top["mask"] == "pLoF"
        assert top["maf"] == "<0.1%"
        assert top["ancestry"] == "All"
        # No integer index and no -log10 may survive into the output.
        assert not any(k.startswith("lp_") for k in top)

    def test_ranks_by_significance(self, pcsk9_gene):
        rows = q.gene_rows(
            pcsk9_gene, ix.phenotypes(), ancestry_idx=0, mask_idx=None,
            maf_idx=None, test="SKAT-O", max_p=None,
        )
        ps = [float(r["p"].lstrip("<")) for r in rows]
        assert ps == sorted(ps)

    def test_max_p_filters(self, pcsk9_gene):
        rows = q.gene_rows(
            pcsk9_gene, ix.phenotypes(), ancestry_idx=0, mask_idx=None,
            maf_idx=None, test="SKAT-O", max_p=1e-50,
        )
        assert rows
        assert all(float(r["p"].lstrip("<")) < 1e-50 for r in rows)

    def test_ci_brackets_the_estimate(self, pcsk9_gene):
        rows = q.gene_rows(
            pcsk9_gene, ix.phenotypes(), ancestry_idx=0, mask_idx=0,
            maf_idx=0, test="SKAT-O", max_p=None,
        )
        lo, hi = (float(x) for x in rows[0]["ci95"].split(" to "))
        assert lo < rows[0]["beta"] < hi


class TestForest:
    def test_returns_every_ancestry_stratum(self, pcsk9_gene):
        res = q.forest(
            pcsk9_gene, ix.phenotypes()[ix.resolve_phenotype("LDLC")],
            pheno_idx=ix.resolve_phenotype("LDLC"), mask_idx=0, maf_idx=0, test="SKAT-O",
        )
        assert [s["ancestry"] for s in res["strata"]] == [
            "All", "EUR", "AFR", "AMR", "EAS", "SAS", "non_EUR"
        ]

    def test_concordance_excludes_pooled_strata(self, pcsk9_gene):
        # 'All' is the comparison baseline and 'non_EUR' pools four of the five
        # superpops, and counting either would double-count the same individuals.
        res = q.forest(
            pcsk9_gene, ix.phenotypes()[ix.resolve_phenotype("LDLC")],
            pheno_idx=ix.resolve_phenotype("LDLC"), mask_idx=0, maf_idx=0, test="SKAT-O",
        )
        assert res["concordance"]["concordant"].endswith("/5")

    def test_concordance_is_flagged_as_derived(self, pcsk9_gene):
        res = q.forest(
            pcsk9_gene, ix.phenotypes()[ix.resolve_phenotype("LDLC")],
            pheno_idx=ix.resolve_phenotype("LDLC"), mask_idx=0, maf_idx=0, test="SKAT-O",
        )
        assert "not a published statistic" in res["concordance"]["basis"]

    def test_unknown_cell_yields_no_strata(self, pcsk9_gene):
        res = q.forest(
            pcsk9_gene, ix.phenotypes()[0], pheno_idx=999, mask_idx=0, maf_idx=0, test="SKAT-O",
        )
        assert res["strata"] == []


class TestCollapse:
    def test_keeps_the_most_significant_row_per_key(self):
        rows = [
            {"ensg": "A", "p": "1e-9"},
            {"ensg": "A", "p": "1e-3"},
            {"ensg": "B", "p": "1e-5"},
        ]
        out = q.collapse_best(rows, ("ensg",))
        assert [r["p"] for r in out] == ["1e-9", "1e-5"]


class TestAggregate:
    def test_ranks_genes_by_how_many_traits_they_hit(self):
        rows = [
            {"ensg": "A", "gene": "AAA", "trait_id": "t1", "trait": "T1", "p": "1e-9"},
            {"ensg": "A", "gene": "AAA", "trait_id": "t2", "trait": "T2", "p": "1e-8"},
            {"ensg": "B", "gene": "BBB", "trait_id": "t1", "trait": "T1", "p": "1e-20"},
        ]
        out = q.aggregate(rows, "gene")
        # Count wins over raw significance: B has the smaller p but hits one trait.
        assert [r["gene"] for r in out] == ["AAA", "BBB"]
        assert out[0]["traits"] == 2
        assert out[0]["strongest_p"] == "1e-9"
        assert out[0]["strongest_trait"] == "T1"
        assert out[0]["trait_list"] == "T1, T2"

    def test_can_rank_traits_by_implicated_genes(self):
        rows = [
            {"ensg": "A", "gene": "AAA", "trait_id": "t1", "trait": "T1", "p": "1e-9"},
            {"ensg": "B", "gene": "BBB", "trait_id": "t1", "trait": "T1", "p": "1e-8"},
            {"ensg": "C", "gene": "CCC", "trait_id": "t2", "trait": "T2", "p": "1e-30"},
        ]
        out = q.aggregate(rows, "trait")
        assert out[0]["trait"] == "T1"
        assert out[0]["genes"] == 2
        assert out[0]["strongest_gene"] == "AAA"


class TestReplicationVerdict:
    """The verdict wording is a scientific claim, so each branch is pinned.

    Conflating "underpowered" with "discordant" is the error that matters here:
    one says the evidence is thin, the other says it points the other way.
    """

    def test_all_strata_significant_is_consistent(self):
        assert q._verdict("5/5", 0.4, "5/5") == "consistent"

    def test_heterogeneity_qualifies_a_consistent_verdict(self):
        assert q._verdict("5/5", 0.001, "5/5") == "consistent but heterogeneous"

    def test_unanimous_direction_with_few_significant_is_underpowered(self):
        assert q._verdict("3/5", 0.4, "5/5") == "same direction in all 5, underpowered in 2"

    def test_genuinely_split_directions_stay_partial(self):
        assert q._verdict("3/5", 0.4, "4/5") == "partial (3/5)"

    def test_nothing_agreeing_is_not_replicated(self):
        assert q._verdict("0/5", 0.4, "2/5") == "not replicated"

    def test_a_missing_count_does_not_invent_a_verdict(self):
        assert q._verdict(None, None, None) == "unknown"
