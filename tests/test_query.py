"""Decoding correctness — the columnar integer encoding is where errors hide."""

import math

import pytest

from brava import index as ix, query as q
from brava.constants import MASK_LABEL, p_from_lp, tier


class TestEffectInterpretation:
    def test_binary_beta_sign_maps_to_risk_or_protection(self):
        assert q.effect_label(0.3, "binary") == "risk+"
        assert q.effect_label(-0.3, "binary") == "protective-"

    def test_quantitative_avoids_risk_language(self):
        # A lower LDL-C is not "protective" in the data's own terms — the trait
        # has no polarity, so the label must not invent one.
        assert q.effect_label(-0.037, "quantitative") == "lower-"
        assert q.effect_label(0.037, "quantitative") == "higher+"

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
        # superpops — counting either would double-count the same individuals.
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
