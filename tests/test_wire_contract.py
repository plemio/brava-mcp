"""Guard on the upstream wire contract.

The published payloads encode every categorical dimension as an integer index
into the canonical lists in `brava.constants`. Upstream documents that contract
as APPEND, NEVER REORDER — but it is their repo, not ours, and a reordering
would not raise: it would silently relabel results, attributing a p-value to the
wrong mask or the wrong ancestry. That is the failure mode worth a network test.

If one of these fails, the fix is to append to `constants` and to open an issue
upstream — never to quietly re-map an index.
"""

import pytest

from brava import client, index as ix
from brava.constants import ANCESTRIES, MAFS, MASKS, TESTS

pytestmark = pytest.mark.network


class TestBundledMetadataStaysInBounds:
    def test_no_unknown_ancestry_appears_in_the_trait_catalogue(self):
        known = set(ANCESTRIES)
        for pheno in ix.phenotypes():
            assert set(pheno["ancestries"]) <= known, pheno["id"]

    @pytest.mark.parametrize("ancestry", ANCESTRIES)
    def test_significant_results_index_decodes_cleanly(self, ancestry):
        payload = ix.all_results(ancestry)
        assert payload["anc"] == ancestry
        assert max(payload["mask_idx"]) < len(MASKS)
        assert max(payload["maf_idx"]) < len(MAFS)
        assert max(payload["test_idx"]) < len(TESTS)
        assert max(payload["gene_idx"]) < len(ix.genes()["ids"])
        assert max(payload["pheno_idx"]) < len(ix.phenotypes())

    def test_gene_index_arrays_stay_aligned(self):
        genes = ix.genes()
        lengths = {len(genes[k]) for k in ("ids", "symbols", "chr", "start", "end")}
        assert len(lengths) == 1


class TestLiveGenePayload:
    async def test_indices_stay_within_our_canonical_lists(self):
        payload = await client.gene_payload("ENSG00000169174")
        assert max(payload["anc"]) < len(ANCESTRIES)
        assert max(payload["mask"]) < len(MASKS)
        assert max(payload["maf"]) < len(MAFS)
        assert max(payload["pheno"]) < len(ix.phenotypes())
        await client.close()

    async def test_every_column_has_one_entry_per_row(self):
        payload = await client.gene_payload("ENSG00000169174")
        cols = ["pheno", "anc", "mask", "maf", "lp_burden", "lp_skat",
                "lp_skato", "lp_het", "beta", "se"]
        assert {len(payload[c]) for c in cols} == {payload["n"]}
        await client.close()
