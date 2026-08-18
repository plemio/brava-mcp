"""End-to-end checks against associations that are not in doubt.

Decoding bugs do not announce themselves: a transposed index yields a
plausible-looking number. These anchor the pipeline on findings whose direction
and rough magnitude are established, so a silent mis-join becomes a test
failure. Mirrors upstream's own pipeline/tests/test_known_biology.py.

Ground truth: https://nikbaya.github.io/brava_browser/
"""

import json

import pytest

import server
from brava import client

pytestmark = pytest.mark.network


def _rows(payload: str) -> list[dict]:
    """The tools emit TOON; re-read it as text for assertions."""
    return payload


class TestPCSK9LowersLDL:
    """PCSK9 loss of function lowers LDL cholesterol. This is the basis of an
    entire drug class; if it comes out positive, the decoding is wrong."""

    async def test_association_is_exome_wide_and_protective(self):
        out = await server.gene_associations("PCSK9", mask="pLoF", maf="<0.1%", limit=1)
        assert "LDL cholesterol" in out
        assert "1.58e-159" in out       # the published meta p-value
        assert "-0.0372" in out          # beta: LOWERS LDL
        assert "exome-wide" in out
        await client.close()

    async def test_direction_holds_in_every_ancestry(self):
        out = await server.gene_phenotype_detail("PCSK9", "LDLC", mask="pLoF", maf="<0.1%")
        assert "5/5" in out              # all five superpopulations concordant
        assert out.count("lowers") == 7  # and all seven strata point the same way
        await client.close()

    async def test_variant_level_reproduces_the_gene_level_direction(self):
        out = await server.gene_variants("PCSK9", "LDLC", max_p=1e-8, limit=5)
        assert "lowers" in out
        assert "concordant" in out
        await client.close()


class TestCanonicalTraitHits:
    async def test_type_2_diabetes_surfaces_the_mody_genes(self):
        # GCK (MODY2) and HNF1A (MODY3) are the textbook rare-variant hits.
        out = await server.phenotype_associations("Type 2 diabetes", limit=5)
        assert "GCK" in out
        assert "HNF1A" in out
        await client.close()

    async def test_lipid_traits_surface_apob(self):
        out = await server.top_associations(category="Cardiovascular", limit=5)
        assert "APOB" in out
        await client.close()


class TestCalibrationControl:
    async def test_synonymous_mask_is_flagged_as_a_control(self):
        # A model must not read a synonymous hit as biology.
        out = await server.gene_associations("PCSK9", mask="synonymous", limit=1)
        assert "CALIBRATION CONTROL" in out
        await client.close()
