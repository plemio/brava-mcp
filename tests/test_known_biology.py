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


class TestTheMaskMafGridIsPinned:
    """Every cell of the (mask, MAF) grid, not just the symmetric one.

    Both dimensions were previously anchored only at index 0 (pLoF, <0.1%),
    which is the one cell where confusing the two columns changes nothing. These
    four cells have four distinct p-values and betas, so any decode that reads
    one dimension from the other's column moves at least one of them.
    """

    @pytest.mark.parametrize(
        "mask,maf,p,beta",
        [
            ("pLoF", "<0.1%", "1.58e-159", "-0.0372"),
            ("pLoF", "<0.01%", "1.12e-76", "-0.0366"),
            ("pLoF | damaging missense", "<0.1%", "1.17e-205", "-0.0267"),
            ("pLoF | damaging missense", "<0.01%", "6.46e-99", "-0.026"),
        ],
    )
    async def test_each_cell_returns_its_own_numbers(self, mask, maf, p, beta):
        out = await server.gene_associations(
            "PCSK9", mask=mask, maf=maf, limit=1, max_p=1e-50
        )
        assert p in out, f"{mask} x {maf} should report p={p}"
        assert beta in out
        await client.close()


class TestCalibrationControl:
    async def test_synonymous_mask_is_flagged_as_a_control(self):
        # A model must not read a synonymous hit as biology.
        out = await server.gene_associations(
            "PCSK9", mask="synonymous", maf="<0.1%", limit=1
        )
        assert "CALIBRATION CONTROL" in out
        await client.close()

    async def test_the_synonymous_row_itself_is_returned(self):
        # The warning above is attached from the mask ARGUMENT, so it would
        # appear even on an empty or wrong result set. Pin the data too: the
        # synonymous control is nowhere near significant, which is the point of
        # a calibration control.
        out = await server.gene_associations(
            "PCSK9", mask="synonymous", maf="<0.1%", limit=1
        )
        assert "0.00324" in out
        assert "ns" in out
        await client.close()
