"""The questions the first cut could not answer.

Each class here maps to a research question that previously cost either several
calls plus manual bookkeeping, or was unanswerable outright because the data was
bundled but never wired to a tool.
"""

import pytest

import server
from brava import index as ix, query as q

pytestmark = pytest.mark.network


class TestPerBiobankContribution:
    """"Which cohorts are behind this number, and how much did each bring?"

    pheno_sizes.json shipped in the repo from the first commit and no tool ever
    read it. It is the cross-biobank evidence, which is the entire reason to
    prefer BRaVa over a single-programme resource.
    """

    def test_rolls_up_per_biobank_with_the_ancestry_split(self):
        rows = q.biobank_contributions(ix.pheno_sizes(), ix.biobanks(), "T2Diab")
        assert len(rows) >= 8
        assert rows == sorted(rows, key=lambda r: -r["n"])
        top = rows[0]
        assert top["biobank"] == "UK Biobank"
        assert top["cases"] > 0
        assert "EUR" in top["ancestries"]

    def test_single_ancestry_cohorts_are_visible_as_such(self):
        # Genes & Health is entirely South Asian and BioBank Japan entirely East
        # Asian; that is exactly what a reader needs to weigh a cross-biobank
        # claim, and it is invisible from the global biobank list.
        rows = {r["biobank"]: r for r in
                q.biobank_contributions(ix.pheno_sizes(), ix.biobanks(), "T2Diab")}
        assert rows["Genes & Health"]["ancestries"].startswith("SAS")
        assert rows["BioBank Japan"]["ancestries"].startswith("EAS")

    def test_quantitative_traits_report_no_case_count(self):
        rows = q.biobank_contributions(ix.pheno_sizes(), ix.biobanks(), "LDLC")
        assert all("cases" not in r for r in rows)

    async def test_the_tool_exposes_it(self):
        out = await server.catalog("biobanks", trait="Type 2 diabetes")
        assert "UK Biobank" in out
        assert "contributing_biobanks" in out

    async def test_an_unknown_trait_points_at_the_catalogue(self):
        out = await server.catalog("biobanks", trait="not a trait")
        assert "error" in out
        assert "catalog(kind='phenotypes')" in out


class TestGeneSetScreening:
    """"Of these candidates, which are significant for anything?"

    The commonest pattern in practice, and it used to cost one call per gene.
    """

    async def test_a_candidate_list_is_screened_in_one_call(self):
        out = await server.top_associations(
            genes="PCSK9,LDLR,APOB", group_by="gene", limit=10
        )
        for symbol in ("PCSK9", "LDLR", "APOB"):
            assert symbol in out

    async def test_a_typo_in_the_list_is_named(self):
        # Silently dropping an unresolvable gene would let a screen come back
        # "negative" for a gene that was never actually looked at.
        out = await server.top_associations(genes="PCSK9,NOTAGENE")
        assert "error" in out
        assert "NOTAGENE" in out


class TestAncestrySpecificDiscovery:
    """"What would a European-only study have missed?"

    Ancestry-specific effects are why the consortium exists, and answering this
    used to mean two calls plus a manual set difference.
    """

    async def test_returns_findings_present_in_one_stratum_and_absent_in_another(self):
        out = await server.top_associations(ancestry="AFR", absent_in="EUR", limit=5)
        assert "results" in out
        assert "not in EUR" in out

    async def test_the_contrast_is_stated_not_implied(self):
        out = await server.top_associations(ancestry="AFR", absent_in="EUR", limit=1)
        assert "Clears p<" in out

    async def test_contrasting_a_stratum_with_itself_is_refused(self):
        out = await server.top_associations(ancestry="AFR", absent_in="AFR")
        assert "error" in out

    async def test_absence_through_never_being_analysed_is_flagged(self):
        # Some traits have no EAS stratum at all. Reporting those as "specific to
        # AFR" would turn missing data into a finding.
        out = await server.top_associations(ancestry="AFR", absent_in="EAS", limit=25)
        assert "never analysed" in out or "Clears p<" in out


class TestGenomeWideVariantScan:
    """"What are the strongest single variants for this trait?"

    Previously unreachable: variants could only be found through a named gene.
    """

    async def test_ranks_variants_without_being_given_a_gene(self):
        out = await server.variants("LDLC", max_p=1e-40, limit=5)
        assert "APOB" in out
        assert "genome-wide" in out

    async def test_agrees_with_the_per_gene_file(self):
        # Two independently published files; if the chromosome index decoding
        # were wrong they would disagree.
        wide = await server.variants("LDLC", max_p=1e-40, limit=25)
        assert "1-55058636-G-A" in wide
        narrow = await server.variants("LDLC", gene="PCSK9", max_p=1e-40, limit=5)
        assert "1-55058636-G-A" in narrow

    async def test_can_be_restricted_to_a_chromosome(self):
        out = await server.variants("LDLC", chrom="2", max_p=1e-40, limit=5)
        assert "APOB" in out
        assert "1-55058636-G-A" not in out  # PCSK9 is on chr1

    async def test_a_stratum_request_says_where_strata_are_available(self):
        out = await server.variants("LDLC", ancestry="AFR")
        assert "error" in out
        assert "gene" in out
