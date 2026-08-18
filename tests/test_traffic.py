"""The traffic commitment made to the upstream author, as executable assertions.

BRaVa's data sits on a personal Cloudflare R2 free tier with hard monthly
ceilings, and upstream's robots.txt asks AI crawlers to stay away. We told them
(nikbaya/brava_browser#1) that search, the catalogue and top-hits queries would
issue NO outbound request at all, and that a given gene would be fetched once
and then never again. A promise nobody checks is not a promise.
"""

import pytest

import server
from brava import client

pytestmark = pytest.mark.network


@pytest.fixture(autouse=True)
def _reset():
    client.reset_counters()
    yield


class TestBundledPathsNeverLeaveTheProcess:
    async def test_search_is_local(self):
        await server.search("PCSK9")
        assert client.outbound_requests() == 0

    async def test_catalogue_is_local(self):
        for kind in ("phenotypes", "biobanks", "vocabulary"):
            await server.catalog(kind)
        assert client.outbound_requests() == 0

    async def test_cross_trait_ranking_is_local(self):
        await server.top_associations(category="Cardiovascular")
        assert client.outbound_requests() == 0

    async def test_trait_top_hits_are_local_at_the_default_threshold(self):
        await server.phenotype_associations("Type 2 diabetes")
        assert client.outbound_requests() == 0


class TestGeneFetchesAreCachedForever:
    async def test_a_cold_gene_costs_exactly_one_request(self):
        # Uses a gene the rest of the suite does not touch, so the disk cache
        # cannot mask a regression here.
        client._memory.clear()
        from brava import index as ix
        ensg = ix.gene_info(ix.resolve_gene("LDLR"))["ensg"]
        client._cache_path(f"{client.DATA_BASE}/gene/{ensg}.json").unlink(missing_ok=True)

        await server.gene_associations("LDLR", limit=1)
        assert client.outbound_requests() == 1

    async def test_a_repeat_lookup_costs_nothing(self):
        await server.gene_associations("PCSK9", limit=1)
        before = client.outbound_requests()
        await server.gene_associations("PCSK9", limit=1)
        await server.gene_phenotype_detail("PCSK9", "LDLC")
        assert client.outbound_requests() == before
