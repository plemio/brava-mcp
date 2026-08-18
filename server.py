"""BRaVa MCP Server.

Gene- and variant-level rare coding-variant association results from the Biobank
Rare Variant Analysis (BRaVa) consortium: ~1.2M individuals, 10 global biobanks,
44 harmonised traits, 7 ancestry strata.

Data source: the public files behind https://nikbaya.github.io/brava_browser/
(MIT, https://github.com/nikbaya/brava_browser). Summary statistics only.

Design notes live in the modules: `index` (bundled metadata + resolution),
`client` (network + traffic discipline), `query` / `variants` (pure decoding).
"""

from __future__ import annotations

import os

import toons
from fastmcp import FastMCP

from brava import client, index as ix, query as q, variants as vq
from brava.constants import (
    ANCESTRIES,
    ANCESTRY_LABEL,
    BROWSER_URL,
    CALIBRATION_MASK_INDEX,
    DEFAULT_ANCESTRY,
    DEFAULT_TEST,
    MAF_LABEL,
    MASK_LABEL,
    MASKS,
    PAPER_URL,
    SIG_GENE_CAUCHY,
    SIG_GENE_MASK_BONFERRONI,
    SIG_SUGGEST,
    SIG_VARIANT,
    TESTS,
    TEST_INDEX,
)

mcp = FastMCP("BRaVa")

DISCLAIMER = "Summary statistics from the BRaVa flagship paper. Not for clinical use."

# Every character of a response crosses the model's context. 25k is the fleet
# budget; an unbounded gene PheWAS at limit=200 with all_tests measured 74k.
CHARACTER_LIMIT = 25_000

# Which parameters would actually shrink an over-budget response, per tool. The
# truncation note names them so the model narrows instead of blindly retrying.
NARROW_GENE = "mask=, maf=, ancestry= or max_p="
NARROW_PHENO = "mask=, maf= or a smaller max_p"

# Read-only, no side effects, and reaching a public dataset outside our control.
READ_ONLY = {"readOnlyHint": True, "idempotentHint": True, "openWorldHint": True}


def _emit(
    payload: dict,
    rows_key: str = "results",
    narrow: str = "",
    offset: int = 0,
) -> str:
    """Serialise, and if the result busts the budget, drop rows and say how to narrow.

    A bare truncation teaches the model nothing and invites a blind retry, so the
    note names the parameters that would actually shrink the answer and the offset
    that continues it. `next_offset` is rewritten to match, otherwise the two
    would disagree and the model would silently skip the dropped rows.

    Row widths vary enough (a 44-trait PheWAS row is not a 3-column catalogue row)
    that a single proportional estimate wastes half the budget, so this binary
    searches for the largest row count that fits.
    """
    out = toons.dumps(payload)
    rows = payload.get(rows_key)
    if len(out) <= CHARACTER_LIMIT or not isinstance(rows, list) or len(rows) <= 1:
        return out

    def render(count: int) -> str:
        trial = {**payload, rows_key: rows[:count]}
        trial["truncated"] = (
            f"Response exceeded the {CHARACTER_LIMIT}-character budget: showing "
            f"{count} of {len(rows)} rows. Continue with offset={offset + count}"
            + (f", or narrow with {narrow}." if narrow else ".")
        )
        if offset + count < offset + len(rows):
            trial["next_offset"] = offset + count
        else:
            trial.pop("next_offset", None)
        return toons.dumps(trial)

    lo, hi, best = 1, len(rows) - 1, render(1)
    while lo <= hi:
        mid = (lo + hi) // 2
        candidate = render(mid)
        if len(candidate) <= CHARACTER_LIMIT:
            best, lo = candidate, mid + 1
        else:
            hi = mid - 1
    return best


def _err(message: str) -> str:
    """Errors are guidance: say what went wrong AND what to do instead."""
    return toons.dumps({"error": message})


def _mask_note(mask_idx: int | None) -> str | None:
    if mask_idx == CALIBRATION_MASK_INDEX:
        return (
            "The 'synonymous' mask is a CALIBRATION CONTROL, not a biological result: "
            "a significant synonymous signal indicates residual test inflation, "
            "not a gene-trait effect."
        )
    return None


# ---------------------------------------------------------------------------

@mcp.tool(annotations=READ_ONLY)
async def search(query: str, limit: int = 10) -> str:
    """Resolve a gene or trait name to the identifiers every other BRaVa tool needs.

    Start here whenever the user names a gene or a disease/trait in free text.
    Accepts a gene symbol ("PCSK9"), an Ensembl id ("ENSG00000169174"), a trait
    abbreviation ("LDLC", "T2Diab"), or a trait name ("LDL cholesterol",
    "Type 2 diabetes"). Matches are ranked exact-id, exact-symbol, prefix, then
    substring.

    Args:
        query: Free text (gene symbol / Ensembl id / trait id / trait name).
        limit: Max results per category (default 10).

    Returns: matching genes (symbol, Ensembl id, GRCh38 position) and matching
    traits (id, name, category, binary/quantitative, sample size).
    """
    limit = max(1, min(limit, 50))
    genes = ix.search_genes(query, limit)
    traits = ix.search_phenotypes(query, limit)
    if not genes and not traits:
        return _err(
            f"No gene or trait matches '{query}'. BRaVa covers 20,033 genes and 44 "
            "harmonised traits. Call catalog(kind='phenotypes') to see the trait list."
        )
    return toons.dumps({"genes": genes, "traits": traits})


@mcp.tool(annotations=READ_ONLY)
async def gene_associations(
    gene: str,
    ancestry: str = DEFAULT_ANCESTRY,
    mask: str | None = None,
    maf: str | None = None,
    test: str = DEFAULT_TEST,
    max_p: float | None = None,
    limit: int = 25,
    offset: int = 0,
    all_tests: bool = False,
) -> str:
    """Phenome-wide association scan for ONE gene: which traits is it associated with?

    Answers "what does rare coding variation in GENE do?" across all 44 BRaVa
    traits, for a chosen ancestry stratum. Results are ranked by p-value.

    Use `gene_phenotype_detail` afterwards to check whether a hit replicates
    across ancestries and biobanks. That cross-biobank view is what BRaVa
    uniquely provides.

    Args:
        gene: Gene symbol or Ensembl id (resolve with `search` if unsure).
        ancestry: All (cross-ancestry meta, default), EUR, AFR, AMR, EAS, SAS,
            non_EUR. Pass "" for every stratum at once.
        mask: Variant annotation mask, one of "pLoF", "damaging missense",
            "other missense", "synonymous", "pLoF | damaging missense",
            "all variants". Omit for all masks.
        maf: Minor-allele-frequency cutoff, "<0.1%" or "<0.01%". Omit for both.
        test: Burden, SKAT, or SKAT-O (default: SKAT-O, the primary omnibus test).
        max_p: Only return associations at or below this p-value.
        limit: Max rows (default 25).
        offset: Skip this many rows, to page through a long result set. The
            response reports total_matching and next_offset.
        all_tests: Also return the other two tests' p-values per row.

    Returns: trait, category, mask, MAF, p-value, significance tier, effect size
    beta with its 95% CI (odds ratio for binary traits), effect direction, and
    the cross-biobank heterogeneity p.
    """
    limit = max(1, min(limit, 200))
    gidx = ix.resolve_gene(gene)
    if gidx is None:
        return _err(f"Unknown gene '{gene}'. Call search('{gene}') to find the right identifier.")

    try:
        anc_idx = ix.resolve_ancestry(ancestry)
        mask_idx = ix.resolve_mask(mask)
        maf_idx = ix.resolve_maf(maf)
        test_name = ix.resolve_test(test, DEFAULT_TEST)
    except ValueError as exc:
        return _err(str(exc))

    info = ix.gene_info(gidx)
    try:
        payload = await client.gene_payload(info["ensg"])
    except client.NotFound:
        return _err(
            f"{info['gene']} ({info['ensg']}) has no BRaVa results. It is in the Ensembl "
            "gene list but was not tested (too few qualifying rare variants)."
        )
    except client.Unavailable as exc:
        return _err(f"BRaVa data is temporarily unreachable: {exc}")

    rows = q.gene_rows(
        payload,
        ix.phenotypes(),
        ancestry_idx=anc_idx,
        mask_idx=mask_idx,
        maf_idx=maf_idx,
        test=test_name,
        max_p=max_p,
        all_tests=all_tests,
    )
    out = {
        "gene": info["gene"],
        "ensg": info["ensg"],
        "position": f"chr{info['chr']}:{info['start']}-{info['end']} (GRCh38)",
        "ancestry": ANCESTRY_LABEL.get(ancestry, "all strata"),
        "test": test_name,
        "total_matching": len(rows),
        "results": rows[offset : offset + limit],
        "note": DISCLAIMER,
    }
    if offset + limit < len(rows):
        out["next_offset"] = offset + limit
    if (mn := _mask_note(mask_idx)) :
        out["warning"] = mn
    return _emit(out, narrow=NARROW_GENE, offset=offset)


@mcp.tool(annotations=READ_ONLY)
async def phenotype_associations(
    phenotype: str,
    ancestry: str = DEFAULT_ANCESTRY,
    mask: str | None = None,
    maf: str | None = None,
    test: str = DEFAULT_TEST,
    max_p: float = SIG_SUGGEST,
    limit: int = 25,
    offset: int = 0,
    detailed: bool = False,
    collapse: bool = True,
) -> str:
    """Top associated genes for ONE trait: which genes carry rare-variant signal?

    The Manhattan view, as a ranked table. By default returns genes clearing the
    suggestive threshold (p < 1e-4) straight from the bundled index, with no
    network access.

    Args:
        phenotype: Trait id or name ("LDLC", "Type 2 diabetes").
        ancestry: All (default), EUR, AFR, AMR, EAS, SAS, non_EUR.
        mask: Variant annotation mask (see `gene_associations`). Omit for all.
        maf: "<0.1%" or "<0.01%". Omit for both.
        test: Burden, SKAT, SKAT-O (default SKAT-O).
        max_p: p-value ceiling (default 1e-4, the suggestive threshold).
        limit: Max rows (default 25).
        offset: Skip this many rows, to page through a long result set.
        detailed: Fetch the full per-trait file to also return standard errors,
            95% CIs and heterogeneity p-values. Required for max_p above 1e-4;
            slower, and the only path that leaves the process.
        collapse: Keep only each gene's most significant mask/MAF combination
            (default true). Set false to see every combination tested.

    Returns: gene symbol, Ensembl id, chromosome, mask, MAF, p-value,
    significance tier and effect size, ranked by p-value.
    """
    limit = max(1, min(limit, 200))
    pidx = ix.resolve_phenotype(phenotype)
    if pidx is None:
        return _err(
            f"Unknown trait '{phenotype}'. Call catalog(kind='phenotypes') for the 44 "
            "traits BRaVa covers, or search('{phenotype}')."
        )
    pheno = ix.phenotypes()[pidx]

    try:
        anc_idx = ix.resolve_ancestry(ancestry)
        mask_idx = ix.resolve_mask(mask)
        maf_idx = ix.resolve_maf(maf)
        test_name = ix.resolve_test(test, DEFAULT_TEST)
    except ValueError as exc:
        return _err(str(exc))
    if anc_idx is None:
        return _err(f"An ancestry is required here. Valid: {', '.join(ANCESTRIES)}")

    anc_name = ANCESTRIES[anc_idx]
    if anc_name not in pheno["ancestries"]:
        return _err(
            f"{pheno['name']} has no {anc_name} stratum. Available: "
            f"{', '.join(pheno['ancestries'])}."
        )

    use_bundle = not detailed and max_p <= SIG_SUGGEST
    if use_bundle:
        rows = q.all_results_rows(
            ix.all_results(anc_name),
            ix.genes(),
            ix.phenotypes(),
            pheno_idx=pidx,
            mask_idx=mask_idx,
            maf_idx=maf_idx,
            test_idx=TEST_INDEX[test_name],
            max_p=max_p,
        )
        source = "bundled significant-results index (p < 1e-4); pass detailed=true for SE, CI and heterogeneity"
    else:
        try:
            payload = await client.phenotype_payload(pheno["id"], anc_name)
        except client.NotFound:
            return _err(f"No BRaVa results file for {pheno['name']} x {anc_name}.")
        except client.Unavailable as exc:
            return _err(f"BRaVa data is temporarily unreachable: {exc}")
        rows = q.phenotype_rows(
            payload,
            ix.genes(),
            pheno["type"],
            mask_idx=mask_idx,
            maf_idx=maf_idx,
            test=test_name,
            max_p=max_p,
        )
        source = "full per-trait results file"

    total = len(rows)
    if collapse:
        rows = q.collapse_best(rows, ("ensg",))

    n_all = (pheno.get("n") or {}).get(anc_name, {})
    out = {
        "trait": pheno["name"],
        "trait_id": pheno["id"],
        "category": pheno["category"],
        "type": pheno["type"],
        "ancestry": ANCESTRY_LABEL[anc_name],
        "sample_size": n_all.get("n"),
        "cases": n_all.get("case"),
        "test": test_name,
        "total_matching": total,
        "distinct_genes": len(rows) if collapse else None,
        "results": rows[offset : offset + limit],
        "source": source
        + ("; showing each gene's most significant mask/MAF only" if collapse else ""),
        "note": DISCLAIMER,
    }
    if offset + limit < len(rows):
        out["next_offset"] = offset + limit
    if (mn := _mask_note(mask_idx)) :
        out["warning"] = mn
    return _emit(out, narrow=NARROW_PHENO, offset=offset)


@mcp.tool(annotations=READ_ONLY)
async def gene_phenotype_detail(
    gene: str,
    phenotype: str,
    mask: str = "pLoF | damaging missense",
    maf: str = "<0.1%",
    test: str = DEFAULT_TEST,
) -> str:
    """Does a gene-trait association replicate across ancestries and biobanks?

    This is BRaVa's distinctive view and the reason to prefer it over a
    single-programme resource. For one gene x trait x mask x MAF cell it returns
    the effect estimate in EVERY ancestry stratum with its 95% CI, the
    cross-cohort heterogeneity p-value, and a concordance tally.

    Use it to qualify any hit found via `gene_associations` or
    `phenotype_associations` before treating it as established.

    Args:
        gene: Gene symbol or Ensembl id.
        phenotype: Trait id or name.
        mask: Variant annotation mask (default "pLoF | damaging missense").
        maf: "<0.1%" (default) or "<0.01%".
        test: Burden, SKAT, SKAT-O (default SKAT-O).

    Returns: one row per ancestry stratum (sample size, p, beta, SE, 95% CI,
    direction), plus heterogeneity and a concordance count. The concordance
    count is DERIVED by this server from upstream's numbers. It is not itself
    a published statistic.
    """
    gidx = ix.resolve_gene(gene)
    if gidx is None:
        return _err(f"Unknown gene '{gene}'. Call search('{gene}').")
    pidx = ix.resolve_phenotype(phenotype)
    if pidx is None:
        return _err(f"Unknown trait '{phenotype}'. Call catalog(kind='phenotypes').")

    try:
        mask_idx = ix.resolve_mask(mask)
        maf_idx = ix.resolve_maf(maf)
        test_name = ix.resolve_test(test, DEFAULT_TEST)
    except ValueError as exc:
        return _err(str(exc))
    if mask_idx is None or maf_idx is None:
        return _err("A specific mask and MAF cutoff are required for this view.")

    info = ix.gene_info(gidx)
    pheno = ix.phenotypes()[pidx]
    try:
        payload = await client.gene_payload(info["ensg"])
    except client.NotFound:
        return _err(f"{info['gene']} ({info['ensg']}) has no BRaVa results.")
    except client.Unavailable as exc:
        return _err(f"BRaVa data is temporarily unreachable: {exc}")

    result = q.forest(
        payload,
        pheno,
        pheno_idx=pidx,
        mask_idx=mask_idx,
        maf_idx=maf_idx,
        test=test_name,
    )
    if not result["strata"]:
        return _err(
            f"No {MASK_LABEL[mask_idx]} / {MAF_LABEL[maf_idx]} result for "
            f"{info['gene']} x {pheno['name']}. Call gene_associations('{info['gene']}') "
            "to see which mask/MAF combinations were tested."
        )

    out = {
        "gene": info["gene"],
        "ensg": info["ensg"],
        "trait": pheno["name"],
        "type": pheno["type"],
        "mask": MASK_LABEL[mask_idx],
        "maf": MAF_LABEL[maf_idx],
        "test": test_name,
        "strata": result["strata"],
        "replication": result["concordance"],
        "note": DISCLAIMER,
    }
    if (mn := _mask_note(mask_idx)) :
        out["warning"] = mn
    return toons.dumps(out)


@mcp.tool(annotations=READ_ONLY)
async def top_associations(
    ancestry: str = DEFAULT_ANCESTRY,
    category: str | None = None,
    trait_type: str | None = None,
    gene: str | None = None,
    mask: str | None = None,
    test: str = DEFAULT_TEST,
    max_p: float = SIG_GENE_CAUCHY,
    limit: int = 25,
    offset: int = 0,
    collapse: bool = True,
    group_by: str | None = None,
) -> str:
    """Strongest rare-variant associations ACROSS traits and genes.

    Answers cross-cutting questions the per-gene and per-trait views cannot:
    "the strongest rare-variant signals in cardiovascular disease", "which
    traits does this gene hit at exome-wide significance", "what are the most
    pleiotropic genes". Served entirely from the bundled index, with no network access.

    Args:
        ancestry: All (default), EUR, AFR, AMR, EAS, SAS, non_EUR.
        category: Trait category filter, e.g. "Cardiovascular", "Lipids",
            "Endocrine/Metabolic" (see catalog(kind='phenotypes')).
        trait_type: "binary" or "quantitative".
        gene: Restrict to one gene. The fast way to ask which traits it hits.
        mask: Variant annotation mask. Omit for all.
        test: Burden, SKAT, SKAT-O (default SKAT-O).
        max_p: p-value ceiling, default 2.5e-6 (the gene-level threshold).
            Cannot exceed 1e-4, the bundled index's inclusion cutoff.
        limit: Max rows (default 25).
        offset: Skip this many rows, to page through a long result set.
        collapse: Return each gene-trait pair once, at its most significant
            mask/MAF (default true). Set false for every combination tested.
        group_by: "gene" ranks genes by how many distinct traits they hit
            (pleiotropy); "trait" ranks traits by how many distinct genes reach
            the threshold. Omit for the flat list of individual associations.

    Returns: gene, trait, category, mask, MAF, test, p-value, tier and effect.
    With group_by, instead returns one row per gene (or trait) with its partner
    count, its strongest partner and p-value, and the partners themselves.
    """
    limit = max(1, min(limit, 200))
    if max_p > SIG_SUGGEST:
        return _err(
            f"max_p above {SIG_SUGGEST} is outside this index, which only holds rows "
            "clearing the suggestive threshold. For weaker signals use "
            "phenotype_associations(detailed=true) on a specific trait."
        )
    try:
        anc_idx = ix.resolve_ancestry(ancestry)
        mask_idx = ix.resolve_mask(mask)
        test_name = ix.resolve_test(test, DEFAULT_TEST)
    except ValueError as exc:
        return _err(str(exc))
    if anc_idx is None:
        return _err(f"An ancestry is required. Valid: {', '.join(ANCESTRIES)}")

    gidx = None
    if gene:
        gidx = ix.resolve_gene(gene)
        if gidx is None:
            return _err(f"Unknown gene '{gene}'. Call search('{gene}').")

    rows = q.all_results_rows(
        ix.all_results(ANCESTRIES[anc_idx]),
        ix.genes(),
        ix.phenotypes(),
        gene_idx=gidx,
        mask_idx=mask_idx,
        test_idx=TEST_INDEX[test_name],
        category=category,
        trait_type=trait_type,
        max_p=max_p,
    )
    if group_by is not None and group_by not in ("gene", "trait"):
        return _err(
            f"Unknown group_by '{group_by}'. Valid: 'gene' (rank genes by how many "
            "traits they hit) or 'trait' (rank traits by how many genes reach the "
            "threshold). Omit it for the flat association list."
        )

    total = len(rows)
    # Aggregation needs one row per pair, so the dedupe is not optional there.
    if collapse or group_by:
        rows = q.collapse_best(rows, ("ensg", "trait_id"))
    pairs = len(rows)
    if group_by:
        rows = q.aggregate(rows, group_by)

    if not rows and category:
        cats = sorted({p["category"] for p in ix.phenotypes()})
        return _err(
            f"No results for category '{category}'. Valid categories: {', '.join(cats)}."
        )
    out = {
        "ancestry": ANCESTRY_LABEL[ANCESTRIES[anc_idx]],
        "test": test_name,
        "p_threshold": max_p,
        "total_matching": total,
        "distinct_pairs": pairs,
        "grouped_by": group_by,
        "results": rows[offset : offset + limit],
        "note": DISCLAIMER,
    }
    if offset + limit < len(rows):
        out["next_offset"] = offset + limit
    return _emit(out, narrow="category=, trait_type=, mask= or a smaller max_p", offset=offset)


@mcp.tool(annotations=READ_ONLY)
async def gene_variants(
    gene: str,
    phenotype: str,
    ancestry: str = DEFAULT_ANCESTRY,
    max_p: float | None = None,
    limit: int = 25,
) -> str:
    """Single-variant results inside a gene, with per-biobank replication.

    Drops from gene-level burden testing down to the individual variants driving
    a signal. Each row carries the per-biobank effect-direction tally, the
    cross-biobank replication evidence no single-programme browser can show,
    plus a gnomAD link for allele frequency, which BRaVa itself does not carry.

    Args:
        gene: Gene symbol or Ensembl id.
        phenotype: Trait id or name.
        ancestry: All (cross-ancestry meta, default) or a specific stratum.
        max_p: p-value ceiling. The variant-level significance threshold is 1.82e-8.
        limit: Max rows (default 25).

    Returns: variant (chr-pos-ref-alt), p-value, beta with 95% CI, effect
    direction, effective sample size, I-squared, heterogeneity p, the ancestries
    it was observed in, the per-biobank direction string, and a gnomAD link.
    """
    limit = max(1, min(limit, 200))
    gidx = ix.resolve_gene(gene)
    if gidx is None:
        return _err(f"Unknown gene '{gene}'. Call search('{gene}').")
    pidx = ix.resolve_phenotype(phenotype)
    if pidx is None:
        return _err(f"Unknown trait '{phenotype}'. Call catalog(kind='phenotypes').")

    try:
        anc_idx = ix.resolve_ancestry(ancestry)
    except ValueError as exc:
        return _err(str(exc))
    if anc_idx is None:
        return _err(f"An ancestry is required. Valid: {', '.join(ANCESTRIES)}")

    info = ix.gene_info(gidx)
    pheno = ix.phenotypes()[pidx]
    split = info["ensg"] in ix.variant_split()
    anc_name = ANCESTRIES[anc_idx]

    try:
        if anc_name == "All":
            payload = await client.gene_variants_payload(
                info["ensg"], pidx if split else None, split
            )
            rows = vq.variant_rows(
                payload, pidx, pheno["type"], max_p=max_p, limit=limit
            )
        else:
            meta = await client.gene_variants_payload(
                info["ensg"], pidx if split else None, split
            )
            payload = await client.gene_variants_anc_payload(
                info["ensg"], pidx if split else None, split
            )
            rows = vq.ancestry_rows(
                payload,
                pidx,
                anc_name,
                pheno["type"],
                max_p=max_p,
                limit=limit,
                chrom=meta.get("chr") or "",
            )
    except client.NotFound:
        return _err(
            f"No variant-level (v2) data for {info['gene']}. Gene-level results are "
            f"still available via gene_associations('{info['gene']}')."
        )
    except client.Unavailable as exc:
        return _err(
            f"Variant-level data is temporarily unreachable: {exc}. Gene-level results "
            f"remain available via gene_associations('{info['gene']}')."
        )
    except (KeyError, TypeError, IndexError) as exc:
        # The v2 format is still evolving upstream; fail legibly, not obscurely.
        return _err(
            f"Variant-level data for {info['gene']} could not be decoded ({exc}). "
            "The upstream v2 format may have changed; gene-level results are unaffected."
        )

    if not rows:
        return toons.dumps(
            {
                "gene": info["gene"],
                "trait": pheno["name"],
                "ancestry": ANCESTRY_LABEL[anc_name],
                "results": [],
                "note": f"No variant-level results for this gene x trait x ancestry. {DISCLAIMER}",
            }
        )

    return _emit(
        {
            "gene": info["gene"],
            "ensg": info["ensg"],
            "trait": pheno["name"],
            "type": pheno["type"],
            "ancestry": ANCESTRY_LABEL[anc_name],
            "variant_significance_threshold": SIG_VARIANT,
            "results": rows,
            "note": f"'biobanks' counts concordant effect directions across contributing "
            f"cohorts; '?' marks a cohort where the variant was absent. {DISCLAIMER}",
        },
        narrow="a smaller max_p",
    )


@mcp.tool(annotations=READ_ONLY)
async def catalog(kind: str = "phenotypes") -> str:
    """What BRaVa covers: traits, contributing biobanks, or the analysis vocabulary.

    Call this before guessing a trait id, or to ground a statement about study
    design (sample sizes, ancestry composition, significance thresholds).

    Args:
        kind: "phenotypes" (the 44 traits with sample sizes per ancestry),
            "biobanks" (the contributing cohorts), or "vocabulary" (masks, MAF
            cutoffs, tests and significance thresholds).

    Returns: the requested catalogue.
    """
    kind = (kind or "phenotypes").strip().lower()

    if kind in ("phenotypes", "phenotype", "traits", "trait"):
        rows = []
        for p in ix.phenotypes():
            n = (p.get("n") or {}).get("All", {})
            rows.append(
                {
                    "trait_id": p["id"],
                    "trait": p["name"],
                    "category": p["category"],
                    "type": p["type"],
                    "n": n.get("n"),
                    "cases": n.get("case"),
                    "ancestries": ",".join(p["ancestries"]),
                    "sex": p.get("sex", ""),
                }
            )
        rows.sort(key=lambda r: (r["category"], r["trait"]))
        return _emit({"traits": rows, "total": len(rows), "note": DISCLAIMER}, "traits")

    if kind in ("biobanks", "biobank", "cohorts"):
        rows = [
            {
                "id": b["id"],
                "name": b["name"],
                "country": b["country"],
                "n": b["sample_size"],
                "sequencing": b["sequencing"],
                "ascertainment": b["ascertainment"],
                "ancestries": ",".join(b["ancestries"]),
            }
            for b in ix.biobanks()
        ]
        rows.sort(key=lambda r: -(r["n"] or 0))
        return _emit({"biobanks": rows, "total": len(rows)}, "biobanks")

    if kind in ("vocabulary", "vocab", "masks", "tests"):
        return toons.dumps(
            {
                "ancestries": [
                    {"id": a, "meaning": ANCESTRY_LABEL[a]} for a in ANCESTRIES
                ],
                "masks": [
                    {"label": MASK_LABEL[i], "raw": MASKS[i]} for i in range(len(MASKS))
                ],
                "maf_cutoffs": MAF_LABEL,
                "tests": [
                    {"name": "SKAT-O", "role": "primary omnibus test; drives significance calls"},
                    {"name": "Burden", "role": "gives the directional effect size (beta) and SE"},
                    {"name": "SKAT", "role": "most sensitive when a gene mixes risk and protective variants"},
                ],
                "significance": {
                    "gene_mask_bonferroni": SIG_GENE_MASK_BONFERRONI,
                    "gene_level_cauchy": SIG_GENE_CAUCHY,
                    "suggestive": SIG_SUGGEST,
                    "variant_level": SIG_VARIANT,
                },
                "caveats": [
                    "The 'synonymous' mask is a calibration control, not a biological result.",
                    "beta > 0 increases risk (binary) or the trait value (quantitative); beta < 0 decreases it.",
                    "beta and SE come from the inverse-variance-weighted Burden meta-analysis, "
                    "even when the reported p-value is SKAT or SKAT-O.",
                    "BRaVa carries no allele frequencies and no common-variant GWAS. "
                    "use gnomAD and the GWAS Catalog for those.",
                ],
                "paper": PAPER_URL,
                "browser": BROWSER_URL,
            }
        )

    return _err(
        f"Unknown catalogue '{kind}'. Valid: 'phenotypes', 'biobanks', 'vocabulary'."
    )


if __name__ == "__main__":
    print("BRaVa MCP Server started", flush=True)
    # Transport switch (same shape as gwas-catalog / ucsc-genome): stdio is the
    # per-worker child the engine spawns; http is the long-lived shared daemon
    # reached via MCP_SERVICES_BASE_URL. Public read-only data, no credentials,
    # so a shared daemon is tenant-safe with no per-request auth.
    transport = os.getenv("MCP_TRANSPORT", "stdio").lower()
    if transport in ("http", "streamable-http"):
        mcp.run(
            transport="http",
            host=os.getenv("MCP_HOST", "127.0.0.1"),
            port=int(os.getenv("MCP_PORT", "3163")),
        )
    else:
        mcp.run()
