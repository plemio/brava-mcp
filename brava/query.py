"""Pure decode / filter / rank / label logic. Stdlib only, no network and no I/O.

The published BRaVa payloads are COLUMNAR: parallel arrays of integer indices
into the canonical lists in `constants`. Nothing in here may leak that encoding
outward: every function returns rows whose values are already labelled and
human-meaningful (a p-value, not a -log10; "pLoF", not 0).

Kept free of aiohttp/fastmcp so the whole decoding surface is unit-testable
against fixtures, which is where the correctness risk actually lives.
"""

from __future__ import annotations

import math
from typing import Any, Iterable

from .constants import (
    ANCESTRIES,
    MAF_LABEL,
    MASK_LABEL,
    SUPERPOPS,
    TEST_LP_KEY,
    p_from_lp,
    tier,
)

# 95% CI multiplier for a normal approximation (what the browser's forest draws).
Z95 = 1.959963984540054

# Keys carrying a p-value in any row this module emits. They are formatted LAST,
# after sorting, because a formatted p is a string and would sort lexically.
P_KEYS = ("p", "p_het", "het_p", "p_burden", "p_skat", "p_skato", "heterogeneity_p")


def fmt_p(p: float | None) -> str | None:
    """p-value as a compact string.

    Serialising 1.58e-159 as a literal decimal costs ~170 characters for one
    cell, the single largest token sink in this server's output, and the reason
    this is not left to the encoder's default float rendering.
    """
    if p is None:
        return None
    if p == 0:
        return "<5e-324"  # upstream floors underflowed p-values; this is real
    if p >= 1e-4:
        return f"{p:.3g}"
    return f"{p:.2e}"


def format_pvalues(rows: Iterable[dict[str, Any]]) -> list[dict[str, Any]]:
    """Format every p-value column in-place. Call only AFTER sorting."""
    out = list(rows)
    for row in out:
        for key in P_KEYS:
            if key in row:
                row[key] = fmt_p(row[key])
    return out


def aggregate(rows: list[dict[str, Any]], by: str) -> list[dict[str, Any]]:
    """Count how many distinct partners each gene (or trait) has, ranked.

    "Which gene is significant for the most traits" and "which trait has the most
    implicated genes" are questions the per-file layout cannot answer at all, and
    that a paginating agent can only answer by pulling every row and tallying by
    hand. Doing it here turns ~16 calls plus manual counting into one call, which
    is the whole point of an outcome-shaped tool.

    Rows must already be deduplicated to one per gene-trait pair and sorted by p
    ascending, so the first row seen for a key is also its strongest.
    """
    key, partner, label = (
        ("ensg", "trait", "traits") if by == "gene" else ("trait_id", "gene", "genes")
    )
    buckets: dict[str, dict[str, Any]] = {}
    for row in rows:
        bucket = buckets.setdefault(
            row[key],
            {
                "gene" if by == "gene" else "trait": row["gene" if by == "gene" else "trait"],
                ("ensg" if by == "gene" else "trait_id"): row[key],
                label: 0,
                "strongest_p": row["p"],
                f"strongest_{partner}": row[partner],
                "partners": [],
            },
        )
        bucket[label] += 1
        bucket["partners"].append(row[partner])

    out = list(buckets.values())
    for bucket in out:
        bucket[partner + "_list"] = ", ".join(bucket.pop("partners")[:12])
    out.sort(key=lambda r: (-r[label], r["strongest_p"]))
    return out


def collapse_best(rows: list[dict[str, Any]], keys: tuple[str, ...]) -> list[dict[str, Any]]:
    """Keep only the most significant row per key tuple.

    Every (gene, trait) pair is tested under 6 masks x 2 MAF cutoffs, so an
    unfiltered ranking returns the same finding a dozen times over and crowds
    out the next real one. Rows MUST already be sorted by p ascending, which
    makes "first seen" identical to "most significant". This is a selection of
    upstream rows, not a new statistic.
    """
    seen: set[tuple] = set()
    out: list[dict[str, Any]] = []
    for row in rows:
        key = tuple(row.get(k) for k in keys)
        if key in seen:
            continue
        seen.add(key)
        out.append(row)
    return out


# ---------------------------------------------------------------------------
# Effect interpretation
# ---------------------------------------------------------------------------

def effect_label(beta: float | None, trait_type: str) -> str | None:
    """Describe an effect direction in terms appropriate to the trait type.

    Binary traits: risk-increasing vs protective. Quantitative traits: the burden
    raises or lowers the trait value, with no good/bad connotation. Mirrors
    upstream app/src/lib/effect.ts.

    Spelled as words rather than a +/- suffix: the sign is already carried by the
    adjacent beta, so a symbol adds nothing and costs a reading step.
    """
    if beta is None or math.isnan(beta) or beta == 0:
        return None
    if trait_type == "binary":
        return "risk-increasing" if beta > 0 else "protective"
    return "raises" if beta > 0 else "lowers"


def ci95(beta: float | None, se: float | None) -> tuple[float, float] | None:
    """Normal-approximation 95% CI around beta."""
    if beta is None or se is None or se <= 0:
        return None
    return (beta - Z95 * se, beta + Z95 * se)


def odds_ratio(beta: float | None, se: float | None) -> dict[str, float] | None:
    """OR = exp(beta) with its CI, for binary traits only.

    Meaningless for quantitative traits, where beta is already on the trait's
    (standardised) scale, so callers must gate on trait type.
    """
    if beta is None:
        return None
    out = {"or": round(math.exp(beta), 4)}
    ci = ci95(beta, se)
    if ci:
        out["or_lo"] = round(math.exp(ci[0]), 4)
        out["or_hi"] = round(math.exp(ci[1]), 4)
    return out


def _sig(x: float | None, digits: int = 3) -> float | None:
    """Round to `digits` significant figures: display precision, not storage."""
    if x is None or x == 0 or math.isnan(x) or math.isinf(x):
        return x
    return round(x, -int(math.floor(math.log10(abs(x)))) + (digits - 1))


# ---------------------------------------------------------------------------
# Row assembly
# ---------------------------------------------------------------------------

def _stat_columns(
    payload: dict,
    i: int,
    test: str,
    trait_type: str,
    all_tests: bool,
) -> dict[str, Any]:
    """The p / beta / effect columns shared by gene- and phenotype-keyed rows."""
    p = p_from_lp(payload[TEST_LP_KEY[test]][i])
    beta = payload["beta"][i]
    se = payload["se"][i]
    row: dict[str, Any] = {"test": test, "p": p, "tier": tier(p)}

    if all_tests:
        for name, key in TEST_LP_KEY.items():
            if name != test:
                row[f"p_{name.lower().replace('-', '')}"] = p_from_lp(payload[key][i])

    row["beta"] = _sig(beta)
    row["se"] = _sig(se)
    ci = ci95(beta, se)
    if ci:
        row["ci95"] = f"{_sig(ci[0])} to {_sig(ci[1])}"
    row["effect"] = effect_label(beta, trait_type)
    if trait_type == "binary":
        orr = odds_ratio(beta, se)
        if orr:
            row["or"] = orr["or"]
            if "or_lo" in orr:
                row["or_ci95"] = f"{orr['or_lo']} to {orr['or_hi']}"
    row["p_het"] = p_from_lp(payload["lp_het"][i])
    return row


def gene_rows(
    payload: dict,
    phenotypes: list[dict],
    *,
    ancestry_idx: int | None,
    mask_idx: int | None,
    maf_idx: int | None,
    test: str,
    max_p: float | None,
    pheno_idx: int | None = None,
    all_tests: bool = False,
) -> list[dict[str, Any]]:
    """Decode `gene/{ENSG}.json` into labelled rows, filtered and ranked by p.

    A row whose selected-test p-value is missing is dropped rather than sorted
    to the end: SAIGE emits null for degenerate strata, and surfacing those as
    "no association" would misread a computational failure as a biological one.
    """
    out: list[dict[str, Any]] = []
    lp_key = TEST_LP_KEY[test]
    for i in range(payload["n"]):
        if pheno_idx is not None and payload["pheno"][i] != pheno_idx:
            continue
        if ancestry_idx is not None and payload["anc"][i] != ancestry_idx:
            continue
        if mask_idx is not None and payload["mask"][i] != mask_idx:
            continue
        if maf_idx is not None and payload["maf"][i] != maf_idx:
            continue
        lp = payload[lp_key][i]
        if lp is None:
            continue
        p = p_from_lp(lp)
        if max_p is not None and p is not None and p > max_p:
            continue

        pheno = phenotypes[payload["pheno"][i]]
        row: dict[str, Any] = {
            "trait": pheno["name"],
            "trait_id": pheno["id"],
            "category": pheno["category"],
            "type": pheno["type"],
            "ancestry": ANCESTRIES[payload["anc"][i]],
            "mask": MASK_LABEL[payload["mask"][i]],
            "maf": MAF_LABEL[payload["maf"][i]],
        }
        row.update(_stat_columns(payload, i, test, pheno["type"], all_tests))
        out.append(row)

    out.sort(key=lambda r: (r["p"] is None, r["p"]))
    return format_pvalues(out)


def phenotype_rows(
    payload: dict,
    genes: dict,
    trait_type: str,
    *,
    mask_idx: int | None,
    maf_idx: int | None,
    test: str,
    max_p: float | None,
    all_tests: bool = False,
) -> list[dict[str, Any]]:
    """Decode `phenotype/{P}.{ANC}.json` into labelled, ranked gene rows.

    The payload is ~2.3 MB of parallel arrays; this is the only thing standing
    between it and the model's context, so filtering happens HERE, never after.
    """
    out: list[dict[str, Any]] = []
    lp_key = TEST_LP_KEY[test]
    for i in range(payload["n"]):
        if mask_idx is not None and payload["mask"][i] != mask_idx:
            continue
        if maf_idx is not None and payload["maf"][i] != maf_idx:
            continue
        lp = payload[lp_key][i]
        if lp is None:
            continue
        p = p_from_lp(lp)
        if max_p is not None and p is not None and p > max_p:
            continue

        g = payload["gene_idx"][i]
        row: dict[str, Any] = {
            "gene": genes["symbols"][g] or genes["ids"][g],
            "ensg": genes["ids"][g],
            "chr": genes["chr"][g],
            "mask": MASK_LABEL[payload["mask"][i]],
            "maf": MAF_LABEL[payload["maf"][i]],
        }
        row.update(_stat_columns(payload, i, test, trait_type, all_tests))
        out.append(row)

    out.sort(key=lambda r: (r["p"] is None, r["p"]))
    return format_pvalues(out)


def all_results_rows(
    payload: dict,
    genes: dict,
    phenotypes: list[dict],
    *,
    pheno_idx: int | None = None,
    gene_idx: int | None = None,
    mask_idx: int | None = None,
    maf_idx: int | None = None,
    test_idx: int | None = None,
    category: str | None = None,
    trait_type: str | None = None,
    max_p: float | None = None,
) -> list[dict[str, Any]]:
    """Decode a bundled `all_results.{ANC}.json` shard.

    This shard holds every row clearing the suggestive cutoff (p < 1e-4) for one
    ancestry, across ALL traits and genes, so it answers cross-trait questions
    the per-file layout cannot, at zero network cost. `beta` is null outside the
    Burden test.
    """
    from .constants import TESTS

    out: list[dict[str, Any]] = []
    for i in range(payload["n"]):
        if pheno_idx is not None and payload["pheno_idx"][i] != pheno_idx:
            continue
        if gene_idx is not None and payload["gene_idx"][i] != gene_idx:
            continue
        if mask_idx is not None and payload["mask_idx"][i] != mask_idx:
            continue
        if maf_idx is not None and payload["maf_idx"][i] != maf_idx:
            continue
        if test_idx is not None and payload["test_idx"][i] != test_idx:
            continue

        pheno = phenotypes[payload["pheno_idx"][i]]
        if category is not None and pheno["category"].lower() != category.lower():
            continue
        if trait_type is not None and pheno["type"] != trait_type:
            continue

        p = p_from_lp(payload["lp"][i])
        if max_p is not None and p is not None and p > max_p:
            continue

        g = payload["gene_idx"][i]
        beta = payload["beta"][i]
        out.append(
            {
                "gene": genes["symbols"][g] or genes["ids"][g],
                "ensg": genes["ids"][g],
                "trait": pheno["name"],
                "trait_id": pheno["id"],
                "category": pheno["category"],
                "ancestry": payload["anc"],
                "mask": MASK_LABEL[payload["mask_idx"][i]],
                "maf": MAF_LABEL[payload["maf_idx"][i]],
                "test": TESTS[payload["test_idx"][i]],
                "p": p,
                "tier": tier(p),
                "beta": _sig(beta),
                "effect": effect_label(beta, pheno["type"]),
            }
        )

    out.sort(key=lambda r: (r["p"] is None, r["p"]))
    return format_pvalues(out)


# ---------------------------------------------------------------------------
# Cross-ancestry forest
# ---------------------------------------------------------------------------

# Strata counted for concordance: the five superpopulations only. 'All' is the
# meta being compared against, and 'non_EUR' pools four of the five. Counting
# either would double-count the same individuals and inflate the tally.
CONCORDANCE_STRATA = SUPERPOPS


def forest(
    payload: dict,
    pheno: dict,
    *,
    pheno_idx: int,
    mask_idx: int,
    maf_idx: int,
    test: str,
) -> dict[str, Any]:
    """Per-ancestry effect estimates for one gene x trait x mask x MAF cell.

    This is the BRaVa-specific view: does the signal replicate across ancestry
    strata, and is it heterogeneous across contributing biobanks? Returns the
    per-stratum numbers plus a concordance tally that is explicitly flagged as
    DERIVED: it is our summary of upstream's numbers, not a published statistic.
    """
    lp_key = TEST_LP_KEY[test]
    by_anc: dict[str, dict[str, Any]] = {}

    for i in range(payload["n"]):
        if payload["pheno"][i] != pheno_idx:
            continue
        if payload["mask"][i] != mask_idx or payload["maf"][i] != maf_idx:
            continue
        anc = ANCESTRIES[payload["anc"][i]]
        beta, se = payload["beta"][i], payload["se"][i]
        p = p_from_lp(payload[lp_key][i])
        row: dict[str, Any] = {
            "ancestry": anc,
            "n": (pheno.get("n") or {}).get(anc, {}).get("n"),
            "p": p,
            "beta": _sig(beta),
            "se": _sig(se),
            "effect": effect_label(beta, pheno["type"]),
        }
        ci = ci95(beta, se)
        if ci:
            row["ci95"] = f"{_sig(ci[0])} to {_sig(ci[1])}"
        if anc == "All":
            row["p_het"] = p_from_lp(payload["lp_het"][i])
        by_anc[anc] = row

    strata = [by_anc[a] for a in ANCESTRIES if a in by_anc]
    meta = by_anc.get("All")

    concordance: dict[str, Any] = {"basis": "derived by this server, not a published statistic"}
    if meta and meta.get("beta"):
        meta_sign = 1 if meta["beta"] > 0 else -1
        counted = [
            r
            for a in CONCORDANCE_STRATA
            if (r := by_anc.get(a)) is not None and r.get("beta") is not None
        ]
        agree = [
            r
            for r in counted
            if (1 if r["beta"] > 0 else -1) == meta_sign
            and r["p"] is not None
            and r["p"] < 0.05
        ]
        same_sign = [r for r in counted if (1 if r["beta"] > 0 else -1) == meta_sign]
        concordance.update(
            {
                "definition": "same effect direction as the meta AND nominal p<0.05, "
                "counted over the 5 superpopulations (All and non_EUR excluded: "
                "they pool the same individuals)",
                "concordant": f"{len(agree)}/{len(counted)}",
                "same_direction": f"{len(same_sign)}/{len(counted)}",
                "heterogeneity_p": meta.get("p_het"),
                "heterogeneity": (
                    "consistent across biobanks"
                    if (meta.get("p_het") or 1) >= 0.05
                    else "heterogeneous across biobanks (p_het < 0.05)"
                ),
            }
        )

    return {"strata": format_pvalues(strata), "concordance": format_pvalues([concordance])[0]}
