"""The four tools, against the real database and the live variant files.

The gene-level surface is now SQL, so most of what used to be a tool parameter
is a clause here. What is still asserted per-tool is what SQL cannot express:
the traps that travel with the schema, the derived replication verdict, and the
output discipline that keeps a 61-million-row table from arriving in the context.
"""

import pytest

import server
from brava import client, db, traps

pytestmark = pytest.mark.network


class TestKnownBiology:
    """Associations that are not in doubt, through the query path.

    A decoding or join error does not announce itself: it yields a plausible
    number. These anchor the pipeline on findings whose direction and magnitude
    are established. Ground truth: https://nikbaya.github.io/brava_browser/
    """

    async def test_pcsk9_lowers_ldl_cholesterol(self):
        # The basis of an entire drug class. A positive beta here means the
        # decoding is wrong, not that the biology changed.
        out = await server.query(
            "SELECT p_skato, beta FROM results WHERE gene='PCSK9' "
            "AND trait_id='LDLC' AND ancestry='All' AND mask='pLoF' AND maf='<0.1%'"
        )
        assert "1.585e-159" in out
        assert "-0.03716" in out

    @pytest.mark.parametrize(
        "mask,maf,p",
        [
            ("pLoF", "<0.1%", "1.585e-159"),
            ("pLoF", "<0.01%", "1.122e-76"),
            ("pLoF | damaging missense", "<0.1%", "1.175e-205"),
            ("pLoF | damaging missense", "<0.01%", "6.457e-99"),
        ],
    )
    async def test_every_cell_of_the_mask_maf_grid_is_distinct(self, mask, maf, p):
        # Anchored off the diagonal on purpose: at (mask=0, maf=0) the two
        # dimensions are indistinguishable, so a join that swaps them passes.
        out = await server.query(
            f"SELECT p_skato FROM results WHERE gene='PCSK9' AND trait_id='LDLC' "
            f"AND ancestry='All' AND mask='{mask}' AND maf='{maf}'"
        )
        assert p in out

    async def test_type_2_diabetes_surfaces_the_mody_genes(self):
        out = await server.query(
            "SELECT gene, min(p_skato) p FROM results WHERE trait_id='T2Diab' "
            "AND ancestry='All' AND mask<>'synonymous' GROUP BY gene ORDER BY p LIMIT 5"
        )
        assert "GCK" in out and "HNF1A" in out

    async def test_the_ancestry_contrast_reproduces_its_published_count(self):
        # Computed independently in evals/resolve_golds.py: 27 distinct pairs.
        out = await server.query(
            "SELECT count(*) n FROM (SELECT DISTINCT a.gene_idx, a.pheno FROM results a "
            "WHERE a.ancestry='AFR' AND a.p_skato < 2.5e-6 AND NOT EXISTS ("
            "SELECT 1 FROM results e WHERE e.ancestry='EUR' AND e.gene_idx=a.gene_idx "
            "AND e.pheno=a.pheno AND e.p_skato < 2.5e-6))"
        )
        assert "27" in out


class TestQueryIsReadOnly:
    @pytest.mark.parametrize(
        "sql",
        ["DROP TABLE assoc", "UPDATE assoc SET beta=0", "ATTACH '/etc/passwd' AS x",
         "SELECT * FROM read_csv('/etc/passwd')", "SELECT 1; SELECT 2", ""],
    )
    async def test_writes_and_file_access_are_refused(self, sql):
        out = await server.query(sql)
        assert "error" in out

    async def test_a_refusal_says_what_to_do_instead(self):
        # A bare rejection invites a blind retry.
        out = await server.query("DELETE FROM assoc")
        assert "schema()" in out or "read statements" in out

    async def test_a_bad_column_points_at_the_schema(self):
        out = await server.query("SELECT nonexistent FROM results LIMIT 1")
        assert "error" in out
        assert "schema()" in out


class TestOutputDiscipline:
    async def test_no_p_value_ever_arrives_as_an_expanded_decimal(self):
        # Four separate builders shipped this bug before the formatting moved to
        # the serialisation boundary, the last being query() itself, where the
        # columns are whatever SQL returned.
        outputs = [
            await server.query("SELECT p_skato, beta FROM results WHERE gene='PCSK9' LIMIT 5"),
            await server.query("SELECT min(p_skato) FROM results"),
            await server.gene_phenotype_detail("PCSK9", "LDLC"),
            await server.schema(),
        ]
        for blob in outputs:
            assert "0.0000000000" not in blob, blob[:200]

    async def test_a_huge_result_set_is_bounded_and_says_so(self):
        out = await server.query("SELECT * FROM results", max_rows=500)
        assert len(out) <= server.CHARACTER_LIMIT
        assert "truncated" in out

    async def test_a_normal_query_is_nowhere_near_the_budget(self):
        out = await server.query(
            "SELECT gene, p_skato FROM results WHERE trait_id='LDLC' "
            "AND ancestry='All' ORDER BY p_skato LIMIT 10"
        )
        assert len(out) < 3000


class TestSchemaCarriesTheTraps:
    """The traps have to arrive with the columns, not in a README nobody fetched."""

    async def test_it_lists_the_shipped_tables_and_the_view_to_use(self):
        out = await server.schema()
        assert "results" in out
        assert "biobank_sizes" in out

    async def test_it_states_where_beta_comes_from(self):
        # A model reading p_skato and beta on the same row will otherwise call it
        # the SKAT-O effect size. There is no such thing.
        out = await server.schema()
        assert "Burden" in out
        assert "no SKAT-O effect size" in out

    async def test_it_flags_the_calibration_control(self):
        assert traps.SYNONYMOUS_IS_A_CONTROL[:40] in await server.schema()

    async def test_it_warns_about_pooled_ancestry_strata(self):
        # Aggregating over every ancestry double-counts: All is the meta and
        # non_EUR pools four of the five superpopulations.
        out = await server.schema()
        assert "double-counts" in out

    async def test_it_explains_that_p_equals_zero_is_the_strongest_result(self):
        out = await server.schema()
        assert "underflow" in out.lower()

    async def test_it_ships_runnable_recipes(self):
        out = await server.schema()
        assert "SELECT" in out
        for question, _ in traps.RECIPES[:2]:
            assert question[:20] in out


class TestReplicationStaysATool:
    """Because the correct query excludes the two pooled strata, and the obvious
    one does not."""

    async def test_a_single_gene_returns_every_stratum(self):
        out = await server.gene_phenotype_detail("PCSK9", "LDLC", mask="pLoF", maf="<0.1%")
        assert out.count("ancestry:") >= 7
        assert "5/5" in out

    async def test_concordance_excludes_the_pooled_strata(self):
        out = await server.gene_phenotype_detail("PCSK9", "LDLC", mask="pLoF", maf="<0.1%")
        assert "/5" in out
        assert "not a published statistic" in out

    async def test_a_list_screens_in_one_call(self):
        out = await server.gene_phenotype_detail("PCSK9,LDLR,APOB", "LDLC")
        assert "screened: 3" in out

    async def test_underpowered_is_not_reported_as_discordant(self):
        out = await server.gene_phenotype_detail("PCSK9,ABCG5", "LDLC")
        assert "same direction in all 5" in out
        assert "underpowered" in out

    async def test_an_unknown_gene_hands_back_a_lookup_query(self):
        out = await server.gene_phenotype_detail("NOTAGENE", "LDLC")
        assert "error" in out
        assert "SELECT" in out


class TestVariantsStayOverHttp:
    async def test_the_genome_wide_scan_needs_no_gene(self):
        out = await server.variants("LDLC", max_p=1e-40, limit=5)
        assert "APOB" in out
        assert "genome-wide" in out

    async def test_within_a_gene_it_reports_per_biobank_concordance(self):
        out = await server.variants("LDLC", gene="PCSK9", max_p=1e-8, limit=3)
        assert "1-55058636-G-A" in out
        assert "concordant" in out

    async def test_an_outage_is_not_reported_as_an_absence(self):
        original = client.variant_overview_payload

        async def unavailable(*args, **kwargs):
            raise client.Unavailable("simulated outage")

        client.variant_overview_payload = unavailable
        try:
            out = await server.variants("LDLC")
        finally:
            client.variant_overview_payload = original
        assert "unreachable" in out
        assert "No genome-wide variant scan published" not in out


class TestTrafficIsNowStructural:
    """Shipping the table changed the traffic story from a promise to a property.

    Gene-level questions used to fetch a file per gene consulted, forever. They
    are local now, so the commitment made upstream (nikbaya/brava_browser#1) is
    no longer something the cache has to keep honouring.
    """

    async def test_the_whole_gene_level_surface_is_local(self):
        client.reset_counters()
        await server.query("SELECT count(*) FROM results WHERE ancestry='All'")
        await server.schema()
        await server.gene_phenotype_detail("PCSK9,LDLR", "LDLC")
        assert client.outbound_requests() == 0

    async def test_only_the_variant_level_leaves_the_process(self):
        client.reset_counters()
        await server.variants("LDLC", max_p=1e-40, limit=1)
        # Cached permanently after the first fetch, so at most one.
        assert client.outbound_requests() <= 1
