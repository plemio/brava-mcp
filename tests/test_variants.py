"""Variant-level (v2) decoding — the densest and least stable wire format."""

from brava import index as ix, variants as vq
from brava.constants import NON_EUR_BIT, SUPERPOPS


class TestAncestryMask:
    def test_decodes_superpop_bits(self):
        assert vq.decode_anc_mask(1 << SUPERPOPS.index("EUR")) == "EUR"
        mask = (1 << SUPERPOPS.index("EUR")) | (1 << SUPERPOPS.index("AFR"))
        assert vq.decode_anc_mask(mask) == "EUR,AFR"

    def test_non_eur_only_is_not_no_data(self):
        # A variant reaching only the pooled non-EUR meta has no superpop bit.
        # Reading that as "no ancestry data" was a real upstream bug (2026-08-17).
        assert vq.decode_anc_mask(NON_EUR_BIT) == "non-EUR (pooled only)"
        assert vq.decode_anc_mask(0) == ""


class TestDirectionSummary:
    def test_counts_the_majority_direction(self):
        assert vq.direction_summary("+++-")["biobanks"] == "3/4 concordant"

    def test_absent_cohorts_do_not_count_as_disagreement(self):
        assert vq.direction_summary("++??")["biobanks"] == "2/2 concordant"

    def test_nothing_to_summarise(self):
        assert vq.direction_summary("????") is None
        assert vq.direction_summary(None) is None


class TestVariantRows:
    def test_joins_slices_onto_the_shared_coordinate_table(self, pcsk9_variants):
        rows = vq.variant_rows(
            pcsk9_variants, ix.resolve_phenotype("LDLC"), "quantitative",
            max_p=None, limit=50,
        )
        assert rows
        chrom, pos, ref, alt = rows[0]["variant"].split("-")
        assert chrom == "1"
        assert int(pos) == rows[0]["pos"]
        assert (ref, alt) == (rows[0]["ref"], rows[0]["alt"])

    def test_links_out_to_gnomad_for_the_frequency_brava_lacks(self, pcsk9_variants):
        rows = vq.variant_rows(
            pcsk9_variants, ix.resolve_phenotype("LDLC"), "quantitative",
            max_p=None, limit=1,
        )
        assert rows[0]["variant"] in rows[0]["gnomad"]

    def test_absent_phenotype_is_an_empty_answer_not_a_crash(self, pcsk9_variants):
        assert vq.variant_rows(pcsk9_variants, 999, "binary", max_p=None, limit=10) == []

    def test_a_missing_column_degrades_to_none(self, pcsk9_variants):
        # Upstream added columns to this format twice in one day; a payload
        # missing one must not raise.
        stripped = {**pcsk9_variants, "by_pheno": {
            k: {kk: vv for kk, vv in v.items() if kk not in ("i2", "ed")}
            for k, v in pcsk9_variants["by_pheno"].items()
        }}
        rows = vq.variant_rows(
            stripped, ix.resolve_phenotype("LDLC"), "quantitative", max_p=None, limit=5
        )
        assert rows and rows[0]["i2"] is None

    def test_a_slice_cannot_outrun_its_coordinate_table(self, pcsk9_variants):
        corrupt = {**pcsk9_variants, "pos": pcsk9_variants["pos"][:2],
                   "ref": pcsk9_variants["ref"][:2], "alt": pcsk9_variants["alt"][:2]}
        rows = vq.variant_rows(
            corrupt, ix.resolve_phenotype("LDLC"), "quantitative", max_p=None, limit=50
        )
        assert all(r["pos"] in corrupt["pos"] for r in rows)
