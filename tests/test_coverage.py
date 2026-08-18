"""The questions the first cut could not answer.

Each class here maps to a research question that previously cost either several
calls plus manual bookkeeping, or was unanswerable outright because the data was
bundled but never wired to a tool.
"""

import pytest

import server
from brava import client, index as ix, query as q

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
        # AFR" would turn missing data into a finding. Asserted on its own, not
        # as an alternative to the always-present contrast line.
        out = await server.top_associations(ancestry="AFR", absent_in="EAS", limit=25)
        assert "never analysed" in out

    async def test_an_outage_is_not_reported_as_an_absence(self):
        # A fetch failure that reads as "this gene has no results" is a claim
        # about biology drawn from a network error.
        original = client.gene_payload

        async def unavailable(*args, **kwargs):
            raise client.Unavailable("simulated outage")

        client.gene_payload = unavailable
        try:
            single = await server.gene_phenotype_detail("PCSK9", "LDLC")
            listed = await server.gene_phenotype_detail("PCSK9,LDLR", "LDLC")
        finally:
            client.gene_payload = original
        for out in (single, listed):
            assert "unreachable" in out
            assert "has no BRaVa results" not in out

    async def test_a_partial_outage_on_non_canonical_input_is_not_double_reported(self):
        # The gene is identified by the resolved index, not by what the caller
        # typed: comparing strings put "pcsk9" in no_result and "PCSK9" in
        # unreachable, so one response said both at once and a model reading the
        # first line concluded biology from a network error.
        original = client.gene_payload

        async def selective(ensg, *args, **kwargs):
            if ensg == "ENSG00000169174":
                raise client.Unavailable("simulated partial outage")
            return await original(ensg, *args, **kwargs)

        client.gene_payload = selective
        try:
            out = await server.gene_phenotype_detail("pcsk9,LDLR", "LDLC")
        finally:
            client.gene_payload = original
        assert "unreachable" in out
        assert "no_result" not in out
        assert "LDLR" in out


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


class TestTheNaiveCallAnswersItsOwnQuestion:
    """"What does PCSK9 do?" used to come back with 4 traits in 25 rows."""

    async def test_a_gene_phewas_reports_distinct_traits(self):
        out = await server.gene_associations("PCSK9", limit=25)
        traits = {l.split(",")[0].strip() for l in out.splitlines() if l.startswith("  ")}
        # Each trait is tested under 6 masks x 2 MAF cutoffs; without collapsing,
        # one trait can occupy half the answer.
        assert len(traits) >= 20
        assert "distinct_traits" in out

    async def test_every_combination_is_still_reachable(self):
        out = await server.gene_associations("PCSK9", collapse=False, limit=12)
        traits = {l.split(",")[0].strip() for l in out.splitlines() if l.startswith("  ")}
        assert len(traits) < 12


class TestBatchReplicationScreen:
    """Qualifying a hit list used to cost one call per gene.

    The tool's own docstring told the model to run it after every hit, so the
    prescribed workflow was a loop: 27 calls and ~44,000 characters to qualify
    LDL cholesterol's hits.
    """

    async def test_a_list_returns_one_verdict_per_gene(self):
        out = await server.gene_phenotype_detail("PCSK9,LDLR,APOB", "LDLC")
        for symbol in ("PCSK9", "LDLR", "APOB"):
            assert symbol in out
        assert "screened: 3" in out
        assert "replication" in out

    async def test_the_screen_is_far_cheaper_than_the_forests(self):
        screen = await server.gene_phenotype_detail("PCSK9,LDLR,APOB,ANGPTL3,ABCG5", "LDLC")
        one = await server.gene_phenotype_detail("PCSK9", "LDLC")
        assert len(screen) < 2 * len(one)

    async def test_a_single_gene_still_gets_the_full_forest(self):
        out = await server.gene_phenotype_detail("PCSK9", "LDLC")
        assert "strata" in out
        assert out.count("ancestry:") >= 7

    async def test_the_verdict_is_labelled_as_derived(self):
        out = await server.gene_phenotype_detail("PCSK9,LDLR", "LDLC")
        assert "not a published statistic" in out

    async def test_underpowered_is_not_reported_as_discordant(self):
        # ABCG5 raises LDL in all five superpopulations but reaches nominal
        # significance in three. Calling that "partial" invites the reader to
        # hear "does not replicate", when the directions are unanimous and two
        # strata are simply small. The distinction is the finding.
        out = await server.gene_phenotype_detail("PCSK9,ABCG5", "LDLC")
        assert "consistent" in out                    # PCSK9: 5/5 significant
        assert "same direction in all 5" in out       # ABCG5: 3/5 significant
        assert "underpowered" in out

    async def test_the_direction_count_is_reported_alongside_the_significant_one(self):
        out = await server.gene_phenotype_detail("ABCG5,PCSK9", "LDLC")
        assert "same_dir" in out

    async def test_an_over_long_list_says_what_to_do(self):
        out = await server.gene_phenotype_detail(",".join(["PCSK9"] * 30), "LDLC")
        assert "error" in out
        assert "top_associations" in out


class TestCandidateScreening:
    async def test_exact_p_for_every_candidate_significant_or_not(self):
        # The screening question is "what is each candidate's p", including the
        # ones that clear nothing. The bundled index only holds p < 1e-4.
        # Deliberately no max_p: the default must not hide a candidate that
        # cleared nothing, which is the whole question a screen asks.
        out = await server.phenotype_associations(
            "LDLC", genes="PCSK9,ACAN", detailed=True, limit=5
        )
        assert "PCSK9" in out and "ACAN" in out
        assert "0.295" in out
        assert "ns" in out

    async def test_a_candidate_with_no_row_at_all_is_named(self):
        # Uses the two Ensembl genes sharing the symbol NOX5, only one of which
        # was tested. Screening by symbol let the untested twin be covered by
        # its namesake and disappear; the previous version of this test passed
        # regardless, because its candidate did have rows.
        out = await server.phenotype_associations(
            "LDLC", genes="ENSG00000290203,ENSG00000255346", detailed=True
        )
        assert "no_result" in out
        assert "ENSG00000290203" in out or "ENSG00000255346" in out

    async def test_the_same_gene_named_twice_is_screened_once(self):
        out = await server.gene_phenotype_detail(
            "PCSK9,ENSG00000169174", "LDLC"
        )
        assert "screened: 2" not in out


class TestProvenance:
    async def test_the_vocabulary_names_the_data_release(self):
        out = await server.catalog("vocabulary")
        assert "data_release" in out
        assert "12 Aug 2026" in out

    async def test_biobanks_report_their_ancestry_composition(self):
        out = await server.catalog("biobanks")
        assert "ancestry_n" in out
        assert "EUR " in out


class TestTruncationIsVisible:
    def test_a_clipped_partner_list_says_so(self):
        rows = [
            {"ensg": "A", "gene": "AAA", "trait_id": f"t{i}", "trait": f"T{i}", "p": "1e-9"}
            for i in range(20)
        ]
        out = q.aggregate(rows, "gene")
        assert out[0]["traits"] == 20
        # A list of 12 next to a count of 20 reads as a contradiction.
        assert "+8 more" in out[0]["trait_list"]
