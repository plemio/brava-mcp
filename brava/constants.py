"""Canonical orderings: a WIRE CONTRACT mirrored from the BRaVa browser ETL.

Integer indices into these lists are the compact keys used throughout the
published JSON, so the order here must match, exactly and in the same order:

    upstream pipeline/common.py   <->   upstream app/src/lib/constants.ts

Upstream documents the contract as APPEND, NEVER REORDER. `tests/test_wire_contract.py`
holds us to that: it re-reads the live metadata and fails if an index we would
decode falls outside these lists. If that test trips, the fix is to append here
(and to open an issue upstream), never to silently re-map.

Source: https://github.com/nikbaya/brava_browser (MIT)
"""

from __future__ import annotations

# --- ancestry strata ---------------------------------------------------------
# 'All' is the cross-ancestry meta-analysis; 'non_EUR' is the pooled non-European
# meta. The suffix is what appears in phenotype filenames ('' for All).
ANCESTRIES: list[str] = ["All", "EUR", "AFR", "AMR", "EAS", "SAS", "non_EUR"]
ANCESTRY_INDEX: dict[str, int] = {a: i for i, a in enumerate(ANCESTRIES)}

ANCESTRY_LABEL: dict[str, str] = {
    "All": "All ancestries (meta-analysis)",
    "EUR": "European",
    "AFR": "African",
    "AMR": "Admixed American",
    "EAS": "East Asian",
    "SAS": "Central & South Asian",
    "non_EUR": "Non-European (meta-analysis)",
}

# The five superpopulations a variant `anc_mask` bitmask can tag, in bit order.
SUPERPOPS: list[str] = ["EUR", "AFR", "AMR", "EAS", "SAS"]
# One bit past the superpops: "also observed in the pooled non-EUR meta".
NON_EUR_BIT = 1 << len(SUPERPOPS)

# --- variant annotation masks ------------------------------------------------
# Raw `Group` strings in canonical index order.
MASKS: list[str] = [
    "pLoF",
    "damaging_missense_or_protein_altering",
    "other_missense_or_protein_altering",
    "synonymous",
    "pLoF;damaging_missense_or_protein_altering",
    "pLoF;damaging_missense_or_protein_altering;other_missense_or_protein_altering;synonymous",
]
MASK_INDEX: dict[str, int] = {m: i for i, m in enumerate(MASKS)}

# Short labels used in output and accepted as input (case-insensitive).
MASK_LABEL: list[str] = [
    "pLoF",
    "damaging missense",
    "other missense",
    "synonymous",
    "pLoF | damaging missense",
    "all variants",
]

# `synonymous` is a CALIBRATION CONTROL, not a biological finding: a significant
# synonymous result signals residual inflation, not a real gene-trait effect.
CALIBRATION_MASK_INDEX = 3

# --- MAF cutoffs -------------------------------------------------------------
MAFS: list[float] = [0.001, 0.0001]
MAF_LABEL: list[str] = ["<0.1%", "<0.01%"]

# --- gene-based tests --------------------------------------------------------
# SKAT-O is the primary omnibus test; Burden is the only one carrying a
# directional effect size (beta/SE); SKAT is most sensitive to mixed directions.
TESTS: list[str] = ["Burden", "SKAT", "SKAT-O"]
TEST_INDEX: dict[str, int] = {t: i for i, t in enumerate(TESTS)}
# Column holding each test's -log10(p) in the gene/phenotype payloads.
TEST_LP_KEY: dict[str, str] = {
    "Burden": "lp_burden",
    "SKAT": "lp_skat",
    "SKAT-O": "lp_skato",
}

# --- significance thresholds (BRaVa flagship paper) --------------------------
SIG_GENE_MASK_BONFERRONI = 1.39e-7  # gene x mask Bonferroni, the strict line
SIG_GENE_CAUCHY = 2.5e-6            # gene-level Cauchy (0.05 / ~20,000 genes)
SIG_SUGGEST = 1e-4                  # suggestive; the all_results inclusion cutoff
SIG_VARIANT = 1.82e-8               # variant-level, 0.05 / 2,746,957

# Smallest positive double: upstream FLOORS underflowed p-values here rather
# than nulling them, so lp == LP_FLOOR means "more significant than we can
# represent", not "missing".
LP_FLOOR = 323.3062153431158

DEFAULT_ANCESTRY = "All"
DEFAULT_TEST = "SKAT-O"
PAPER_URL = "https://www.medrxiv.org/content/10.64898/2026.05.21.26353759v1.full"
BROWSER_URL = "https://nikbaya.github.io/brava_browser/"


def tier(p: float | None) -> str:
    """Significance tier for a p-value, using the paper's own thresholds."""
    if p is None:
        return "n/a"
    if p < SIG_GENE_MASK_BONFERRONI:
        return "exome-wide"
    if p < SIG_GENE_CAUCHY:
        return "gene-level"
    if p < SIG_SUGGEST:
        return "suggestive"
    return "ns"


def p_from_lp(lp: float | None) -> float | None:
    """-log10(p) -> p. Returns 0.0 at the floor (a real underflow, not missing)."""
    if lp is None:
        return None
    if lp >= LP_FLOOR:
        return 0.0
    return 10.0**-lp
