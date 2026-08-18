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
P_KEYS = (
    "p",
    "p_het",
    "het_p",
    "p_burden",
    "p_skat",
    "p_skato",
    "heterogeneity_p",
    "strongest_p",  # produced by aggregate(), which runs before formatting
)


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

    Rows must arrive UNFORMATTED (numeric `p`) and already deduplicated to one row
    per gene-trait pair, sorted by p ascending, so the first row seen for a key is
    also its strongest. Formatted p-values are strings, and sorting those compares
    them lexically: "1.00e-08" would outrank "9.99e-300", and "<5e-324" (the most
    significant value there is) would sort last of all. `format_pvalues` is called
    on the RESULT of this function, never before it.
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
        names = bucket.pop("partners")
        shown = names[:12]
        # Say when the list is clipped: an incomplete list sitting next to a
        # complete count reads as a contradiction, or worse, as the whole set.
        bucket[partner + "_list"] = ", ".join(shown) + (
            f", +{len(names) - len(shown)} more" if len(names) > len(shown) else ""
        )
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
    # Emitted for every row, empty where SE is missing or zero. SAIGE produces
    # degenerate rows with se=0 at the p-value floor, which sort FIRST, so a
    # conditional key here is not a rare edge: APOB pLoF <0.01% has four of them
    # and the response fell out of TOON's table form, 4,740 characters to 8,895.
    ci = ci95(beta, se)
    row["ci95"] = f"{_sig(ci[0])} to {_sig(ci[1])}" if ci else ""
    row["effect"] = effect_label(beta, trait_type)
    # Emitted for every row, empty on quantitative traits where an odds ratio is
    # meaningless. Keys must not vary between rows: TOON encodes a uniform list
    # as one compact table and a ragged one as field-per-line, and a single
    # optional key doubled the size of a mixed-trait response.
    orr = odds_ratio(beta, se) if trait_type == "binary" else None
    row["or"] = orr["or"] if orr else ""
    row["or_ci95"] = (
        f"{orr['or_lo']} to {orr['or_hi']}" if orr and "or_lo" in orr else ""
    )
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
    gene_idxs: set[int] | None = None,
) -> list[dict[str, Any]]:
    """Decode `phenotype/{P}.{ANC}.json` into labelled, ranked gene rows.

    The payload is ~2.3 MB of parallel arrays; this is the only thing standing
    between it and the model's context, so filtering happens HERE, never after.
    """
    out: list[dict[str, Any]] = []
    lp_key = TEST_LP_KEY[test]
    for i in range(payload["n"]):
        if gene_idxs is not None and payload["gene_idx"][i] not in gene_idxs:
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
    gene_idxs: set[int] | None = None,
    mask_idx: int | None = None,
    maf_idx: int | None = None,
    test_idx: int | None = None,
    category: str | None = None,
    trait_type: str | None = None,
    max_p: float | None = None,
    format_p: bool = True,
) -> list[dict[str, Any]]:
    """Decode a bundled `all_results.{ANC}.json` shard.

    This shard holds every row clearing the suggestive cutoff (p < 1e-4) for one
    ancestry, across ALL traits and genes, so it answers cross-trait questions
    the per-file layout cannot, at zero network cost.

    `beta` and `se` always come from the inverse-variance-weighted Burden
    meta-analysis, including on SKAT and SKAT-O rows: all 4,348 SKAT-O rows in the
    All shard carry one. (Upstream's types.ts says beta is null outside Burden;
    the published data disagrees, so the data wins.)

    Set format_p=False to keep p numeric, which callers that sort or aggregate
    downstream must do.
    """
    from .constants import TESTS

    out: list[dict[str, Any]] = []
    for i in range(payload["n"]):
        if pheno_idx is not None and payload["pheno_idx"][i] != pheno_idx:
            continue
        if gene_idx is not None and payload["gene_idx"][i] != gene_idx:
            continue
        if gene_idxs is not None and payload["gene_idx"][i] not in gene_idxs:
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
    return format_pvalues(out) if format_p else out


def replication_summary(forests: list[tuple[str, dict]]) -> list[dict[str, Any]]:
    """Compact replication verdict per gene, for screening a hit list.

    The full forest is seven rows of estimates per gene, which is right for one
    gene and unusable for twenty-seven: qualifying every LDL-C hit that way costs
    27 calls and ~44k characters. This keeps only what the screening question
    turns on, so a whole hit list fits in one table, and the caller drops back to
    the full forest for whichever rows deserve a closer look.
    """
    out: list[dict[str, Any]] = []
    for gene, forest_result in forests:
        strata = {row["ancestry"]: row for row in forest_result["strata"]}
        meta = strata.get("All", {})
        conc = forest_result["concordance"]
        agree = conc.get("concordant", "")
        het_p = conc.get("heterogeneity_p")
        out.append(
            {
                "gene": gene,
                "p": meta.get("p"),
                "tier": tier_of(meta.get("p")),
                "beta": meta.get("beta"),
                "effect": meta.get("effect"),
                "concordant": agree,
                "same_dir": conc.get("same_direction", ""),
                "p_het": het_p,
                "replication": _verdict(agree, het_p, conc.get("same_direction")),
                "strata_tested": len(strata),
            }
        )
    return out


def tier_of(p: str | float | None) -> str:
    """Tier from an already-formatted p-value string (or a float)."""
    if p is None:
        return "n/a"
    if isinstance(p, str):
        p = 0.0 if p.startswith("<") else float(p)
    return tier(p)


def _verdict(
    concordant: str,
    het_p: str | float | None,
    same_direction: str | None = None,
) -> str:
    """One word for "did this hold up", derived and labelled as such.

    Separates underpowered from discordant, because collapsing them is a real
    scientific error and not a rounding of nuance: ANGPTL3 x LDL-C lowers LDL in
    all five superpopulations, but AFR (p=0.11) and SAS miss nominal
    significance on sample size alone. Reporting that as "partial" invites the
    reader to hear "does not replicate", when the directions are unanimous.

    Deliberately coarse otherwise. The inputs are upstream's; the grouping is
    ours, which is why the caller states it is derived rather than published.
    """
    def _split(value: str | None) -> tuple[int, int] | None:
        try:
            a, b = (int(x) for x in str(value).split("/"))
            return a, b
        except (ValueError, AttributeError):
            return None

    parsed = _split(concordant)
    if parsed is None:
        return "unknown"
    agree, total = parsed
    if total == 0:
        # No stratum carries an estimate: nothing was compared, so neither
        # "replicated" nor "not replicated" is a statement the data supports.
        return "inconclusive (no stratum estimates)"
    same = _split(same_direction)

    het = None
    if het_p is not None:
        het = 0.0 if str(het_p).startswith("<") else float(het_p)
    heterogeneous = het is not None and het < 0.05

    if total and agree == total:
        return "consistent but heterogeneous" if heterogeneous else "consistent"
    if same and same[0] == same[1]:
        base = f"same direction in all {same[1]}, underpowered in {same[1] - agree}"
        return base + (", heterogeneous" if heterogeneous else "")
    if agree == 0:
        # Nothing reaches nominal significance anywhere. "Not replicated" would
        # assert a negative the data cannot carry; it was never tested with power.
        return "inconclusive (no stratum reaches nominal significance)"
    return f"partial ({concordant})" + (", heterogeneous" if heterogeneous else "")


def significant_pairs(payload: dict, max_p: float, test_idx: int) -> set[tuple[int, int]]:
    """(gene, trait) pairs clearing `max_p` in one ancestry shard.

    The building block for "significant here but not there". The shards only hold
    rows past the suggestive cutoff, so absence from this set has THREE possible
    causes, and only the caller can tell them apart:

      1. tested in that stratum, did not clear `max_p`;
      2. the stratum was never analysed for that trait (the caller checks the
         trait's ancestry list);
      3. the gene had too few qualifying carriers to be tested at all in that
         stratum, which is common outside Europe and indistinguishable from (1)
         from this side.

    Cause 3 is why a contrast is never on its own evidence of an ancestry-specific
    effect: the comparison strata differ in size by an order of magnitude, so a
    difference in significance is not a significant difference.
    """
    out: set[tuple[int, int]] = set()
    for i in range(payload["n"]):
        if payload["test_idx"][i] != test_idx:
            continue
        if p_from_lp(payload["lp"][i]) > max_p:
            continue
        out.add((payload["gene_idx"][i], payload["pheno_idx"][i]))
    return out


def biobank_contributions(
    sizes: dict,
    biobanks: list[dict],
    trait_id: str,
) -> list[dict[str, Any]]:
    """Who actually contributed to one trait's analysis, and how much.

    BRaVa's whole claim is that a finding held across ten independent biobanks,
    so "which cohorts are behind this number" is not trivia: it is the question
    that separates a replicated result from one carried by a single cohort.
    Rolls the per-(trait, superpopulation, biobank) sizes up per biobank, keeping
    the ancestry breakdown inline rather than as 50 near-empty rows.
    """
    by_pop = sizes.get(trait_id) or {}
    names = {b["id"]: b for b in biobanks}
    rolled: dict[str, dict[str, Any]] = {}

    for pop, entries in by_pop.items():
        for entry in entries:
            row = rolled.setdefault(
                entry["id"],
                {
                    "biobank": names.get(entry["id"], {}).get("name", entry["id"]),
                    "country": names.get(entry["id"], {}).get("country", ""),
                    "n": 0,
                    "cases": 0,
                    "_pops": [],
                },
            )
            row["n"] += entry.get("n") or 0
            row["cases"] += entry.get("case") or 0
            row["_pops"].append(f"{pop} {entry.get('n') or 0:,}")

    out = list(rolled.values())
    for row in out:
        row["ancestries"] = "; ".join(sorted(row.pop("_pops")))
        if not row["cases"]:
            # Quantitative trait: a case count of zero is absence, not a measurement.
            row.pop("cases")
    out.sort(key=lambda r: -r["n"])
    return out


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
        # Unconditional, like every other optional column: a key that appears
        # on some rows and not others drops TOON out of its table encoding.
        ci = ci95(beta, se)
        row["ci95"] = f"{_sig(ci[0])} to {_sig(ci[1])}" if ci else ""
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
