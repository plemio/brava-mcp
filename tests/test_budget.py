"""Output budget and pagination.

Every character of a response crosses the model's context. An unbounded gene
PheWAS at limit=200 with all_tests measured 74k characters, roughly 18k tokens
for a single call, so the budget is enforced rather than trusted.
"""

import pytest
import toons

import server
from brava import client, query as q

pytestmark = pytest.mark.network


class TestCharacterBudget:
    async def test_a_deliberately_huge_request_is_bounded(self):
        out = await server.gene_associations(
            "PCSK9", ancestry="", limit=200, all_tests=True
        )
        assert len(out) <= server.CHARACTER_LIMIT

    async def test_truncation_names_the_parameters_that_would_narrow_it(self):
        # A bare cut teaches nothing and invites a blind retry.
        out = await server.gene_associations(
            "PCSK9", ancestry="", limit=200, all_tests=True
        )
        assert "truncated" in out
        assert "mask=" in out
        assert "offset=" in out

    async def test_the_continuation_offset_agrees_with_the_truncation(self):
        # If next_offset kept its pre-truncation value the model would page past
        # the rows that were dropped, losing them silently.
        out = await server.gene_associations(
            "PCSK9", ancestry="", limit=200, offset=50, all_tests=True
        )
        shown = int(out.split("showing ")[1].split(" of ")[0])
        stated = int(out.split("Continue with offset=")[1].split(",")[0])
        next_offset = int(
            [l for l in out.splitlines() if l.startswith("next_offset:")][0].split(":")[1]
        )
        assert stated == 50 + shown == next_offset

    async def test_a_normal_call_is_nowhere_near_the_budget(self):
        out = await server.gene_associations("PCSK9")
        assert len(out) < 6000
        assert "truncated" not in out


class TestPagination:
    async def test_offset_advances_the_window(self):
        first = await server.top_associations(limit=3)
        second = await server.top_associations(limit=3, offset=3)
        assert first != second

    async def test_next_offset_appears_only_while_rows_remain(self):
        many = await server.top_associations(limit=3)
        assert "next_offset" in many
        # 44 traits, so a limit far past the total must not advertise another page.
        few = await server.catalog("phenotypes")
        assert "next_offset" not in few

    async def test_paging_past_the_end_is_empty_not_an_error(self):
        out = await server.gene_associations("PCSK9", offset=100_000)
        assert "error" not in out
        assert "results[0]" in out or "results:" in out


class TestReadOnlyAnnotations:
    async def test_every_tool_declares_itself_read_only(self):
        # Nothing here writes anything, and all of it reaches a public dataset
        # outside our control.
        tools = await server.mcp.list_tools()
        assert len(tools) == 7
        for t in tools:
            assert t.annotations.readOnlyHint is True, t.name
            assert t.annotations.openWorldHint is True, t.name


class TestAggregationPath:
    async def test_pleiotropy_is_one_call_not_sixteen(self):
        # Without group_by a model must page 388 rows at limit=25 and tally by
        # hand. That is the "correct in 25 calls" failure the design avoids.
        out = await server.top_associations(max_p=1.39e-7, group_by="gene", limit=3)
        assert "PKD1" in out
        assert "traits" in out
        assert client.outbound_requests() == 0

    async def test_an_unknown_grouping_names_the_valid_ones(self):
        out = await server.top_associations(group_by="chromosome")
        assert "error" in out
        assert "'gene'" in out and "'trait'" in out


class TestNoRawPValueEverEscapes:
    """A p-value serialised as a literal decimal costs ~170 characters.

    Three separate code paths have shipped this bug now, each time because a new
    row builder forgot the formatting pass. This asserts the property across
    every tool at once, so the next builder is covered whether or not whoever
    writes it remembers.
    """

    async def test_no_tool_emits_an_expanded_decimal(self):
        outputs = [
            await server.gene_associations("PCSK9", limit=3),
            await server.phenotype_associations("LDLC", limit=3),
            await server.gene_phenotype_detail("PCSK9", "LDLC"),
            await server.top_associations(limit=3),
            await server.top_associations(group_by="gene", limit=3),
            await server.variants("LDLC", max_p=1e-40, limit=3),
            await server.variants("LDLC", gene="PCSK9", limit=3),
            await server.catalog("vocabulary"),
        ]
        for blob in outputs:
            assert "0.0000000000" not in blob, blob[:200]


class TestRowsKeepAUniformShape:
    """Optional keys are a token bomb, not a cosmetic detail.

    TOON encodes a uniform list of rows as one table with a single header, and a
    ragged one field-per-line. A mixed binary/quantitative response where only
    the binary rows carried `or` doubled in size, which surfaced the moment
    collapsing made a single response span many traits.
    """

    async def test_a_mixed_trait_response_stays_tabular(self):
        out = await server.gene_associations("PCSK9", limit=25)
        # The table header appears once; a ragged encoding has no header at all.
        assert out.count("results[") == 1
        assert "{trait,trait_id" in out

    async def test_mixing_binary_and_quantitative_costs_nothing(self):
        mixed = await server.gene_associations("PCSK9", limit=25)
        single = await server.gene_associations("PCSK9", collapse=False, limit=25)
        assert len(mixed) < 1.5 * len(single)

    async def test_a_degenerate_row_does_not_break_the_table(self):
        # SAIGE emits rows with se=0 at the p-value floor, and those sort FIRST.
        # APOB pLoF <0.01% has four; while ci95 was conditional this response
        # fell out of the table form, 4,740 characters to 8,895.
        out = await server.gene_associations("APOB", mask="pLoF", maf="<0.01%")
        assert "{trait,trait_id" in out
        assert len(out) < 6000

    def test_every_optional_column_is_present_on_a_degenerate_row(self):
        from brava import index as ix

        rows = q.gene_rows(
            {
                "n": 1, "pheno": [0], "anc": [0], "mask": [0], "maf": [0],
                "lp_burden": [400.0], "lp_skat": [400.0], "lp_skato": [400.0],
                "lp_het": [None], "beta": [-0.5], "se": [0.0],
            },
            ix.phenotypes(),
            ancestry_idx=0, mask_idx=0, maf_idx=0, test="SKAT-O", max_p=None,
        )
        # Every key a healthy row carries must be here, empty rather than absent.
        for key in ("ci95", "or", "or_ci95"):
            assert key in rows[0], key

    def test_odds_ratio_columns_exist_even_when_empty(self):
        from brava import index as ix

        rows = q.gene_rows(
            {
                "n": 1, "pheno": [0], "anc": [0], "mask": [0], "maf": [0],
                "lp_burden": [5.0], "lp_skat": [5.0], "lp_skato": [5.0],
                "lp_het": [1.0], "beta": [0.1], "se": [0.01],
            },
            ix.phenotypes(),
            ancestry_idx=0, mask_idx=0, maf_idx=0, test="SKAT-O", max_p=None,
        )
        assert "or" in rows[0] and "or_ci95" in rows[0]
