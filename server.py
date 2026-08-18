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

import asyncio
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

# beta/SE come from the Burden meta-analysis even on rows whose reported p-value
# is SKAT or SKAT-O. Stated per response rather than only in catalog(vocabulary):
# a client that never calls the catalogue would otherwise read the beta as
# belonging to the test named in the same row.
BETA_NOTE = "beta/SE: inverse-variance-weighted Burden meta-analysis, whichever test the p-value comes from."

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


def _ambiguity_note(query: str, chosen: int) -> str | None:
    """Warn when a gene symbol maps to more than one Ensembl gene."""
    others = ix.ambiguous_alternatives(query, chosen)
    if not others:
        return None
    return (
        f"'{query}' maps to more than one Ensembl gene. These results are for "
        f"{ix.gene_info(chosen)['ensg']}; also matching: {', '.join(others)}. "
        "Pass an Ensembl id to select a specific one."
    )


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
    collapse: bool = True,
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
        collapse: Report each trait once, at the mask/MAF where it is most
            significant (default true). Each trait is tested under 6 masks x 2
            MAF cutoffs, so without this a 25-row answer covers only a handful
            of traits. Set false to see every combination tested.

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
    total = len(rows)
    if collapse:
        # A trait repeated across 12 mask/MAF cells crowds out the next trait,
        # so the flagship "what does this gene do" call answered with 4 distinct
        # traits out of 25 rows.
        rows = q.collapse_best(rows, ("trait_id", "ancestry"))
    out = {
        "gene": info["gene"],
        "ensg": info["ensg"],
        "position": f"chr{info['chr']}:{info['start']}-{info['end']} (GRCh38)",
        "ancestry": ANCESTRY_LABEL.get(ancestry, "all strata"),
        "test": test_name,
        "total_matching": total,
        # trait x stratum when no ancestry filter is applied: there are only 44
        # traits, so calling 280 rows "distinct traits" is a false label.
        ("distinct_traits" if anc_idx is not None else "distinct_trait_strata"): (
            len(rows) if collapse else None
        ),
        "results": rows[offset : offset + limit],
        "note": f"{BETA_NOTE} {DISCLAIMER}",
    }
    if offset + limit < len(rows):
        out["next_offset"] = offset + limit
    notes = [n for n in (_ambiguity_note(gene, gidx), _mask_note(mask_idx)) if n]
    if notes:
        out["warning"] = " ".join(notes)
    return _emit(out, narrow=NARROW_GENE, offset=offset)


@mcp.tool(annotations=READ_ONLY)
async def phenotype_associations(
    phenotype: str,
    ancestry: str = DEFAULT_ANCESTRY,
    mask: str | None = None,
    maf: str | None = None,
    test: str = DEFAULT_TEST,
    max_p: float | None = None,
    limit: int = 25,
    offset: int = 0,
    detailed: bool = False,
    collapse: bool = True,
    genes: str | None = None,
    all_tests: bool = False,
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
        max_p: p-value ceiling. Defaults to 1e-4 (the suggestive threshold) when
            ranking a whole trait, and to NO ceiling when `genes` is given: a
            candidate screen that hides the candidates which cleared nothing
            answers the opposite of the question asked.
        limit: Max rows (default 25).
        offset: Skip this many rows, to page through a long result set.
        detailed: Fetch the full per-trait file to also return standard errors,
            95% CIs and heterogeneity p-values. Required for max_p above 1e-4;
            slower, and the only path that leaves the process.
        collapse: Keep only each gene's most significant mask/MAF combination
            (default true). Set false to see every combination tested.
        genes: Restrict to a comma-separated candidate list ("PCSK9,LDLR,APOB").
            With detailed=true this reports each candidate's exact p-value even
            where it clears no threshold, which is the screening question, in one
            call instead of one per gene.
        all_tests: Also return the other two tests' p-values per row.

    Returns: gene symbol, Ensembl id, chromosome, mask, MAF, p-value,
    significance tier and effect size, ranked by p-value.
    """
    limit = max(1, min(limit, 200))
    pidx = ix.resolve_phenotype(phenotype)
    if pidx is None:
        return _err(
            f"Unknown trait '{phenotype}'. Call catalog(kind='phenotypes') for the 44 "
            f"traits BRaVa covers, or search('{phenotype}')."
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

    # Resolved here rather than in the signature: the right default depends on
    # whether this is a ranking or a screen, and a plain default cannot tell
    # "not passed" from "passed the same value".
    screening = bool(genes)
    if max_p is None:
        max_p = 1.0 if screening else SIG_SUGGEST

    gene_idxs = None
    wanted: list[str] = []
    resolved: dict[str, int | None] = {}
    if genes:
        wanted = [g.strip() for g in genes.split(",") if g.strip()]
        resolved = {g: ix.resolve_gene(g) for g in wanted}
        unknown = [g for g, i in resolved.items() if i is None]
        if unknown:
            return _err(
                f"Unknown gene(s): {', '.join(unknown)}. Call search() on each, or "
                "pass Ensembl ids. Screening a list is only meaningful if every "
                "name resolved."
            )
        gene_idxs = set(resolved.values())

    # A candidate screen wants each gene's p-value whether or not it is
    # significant, and the bundled index only holds rows past 1e-4.
    use_bundle = not detailed and max_p <= SIG_SUGGEST and gene_idxs is None
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
            all_tests=all_tests,
            gene_idxs=gene_idxs,
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
        "note": f"{BETA_NOTE} {DISCLAIMER}",
    }
    if offset + limit < len(rows):
        out["next_offset"] = offset + limit
    if screening:
        shown = {r["gene"] for r in rows} | {r["ensg"] for r in rows}
        absent = [
            g for g in wanted
            if g not in shown and ix.gene_info(resolved[g])["gene"] not in shown
        ]
        if absent:
            # Never let a candidate vanish from a screen: silence reads as "not
            # associated", the exact false negative this is meant to prevent.
            out["no_result"] = (
                f"{', '.join(absent)}: no row under the requested mask/MAF filters. "
                "That is not the same as tested and null."
            )
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
    `phenotype_associations` before treating it as established. Pass a
    comma-separated list to screen a whole hit list at once: qualifying every
    LDL-C hit one at a time costs 27 calls and ~44,000 characters, so the list
    form returns a compact verdict per gene and you come back here for whichever
    rows deserve the full forest.

    Args:
        gene: Gene symbol or Ensembl id, or a comma-separated list of them
            ("PCSK9,LDLR,APOB") to screen several at once.
        phenotype: Trait id or name.
        mask: Variant annotation mask (default "pLoF | damaging missense").
        maf: "<0.1%" (default) or "<0.01%".
        test: Burden, SKAT, SKAT-O (default SKAT-O).

    Returns: one row per ancestry stratum (sample size, p, beta, SE, 95% CI,
    direction), plus heterogeneity and a concordance count. The concordance
    count is DERIVED by this server from upstream's numbers. It is not itself
    a published statistic.
    """
    wanted = [g.strip() for g in (gene or "").split(",") if g.strip()]
    if not wanted:
        return _err("Pass a gene symbol or Ensembl id, or a comma-separated list.")
    if len(wanted) > 25:
        return _err(
            f"{len(wanted)} genes is more than this screens at once (max 25). "
            "Split the list, or narrow it with top_associations first."
        )
    resolved = {g: ix.resolve_gene(g) for g in wanted}
    unknown = [g for g, i in resolved.items() if i is None]
    if unknown:
        return _err(f"Unknown gene(s): {', '.join(unknown)}. Call search() on each.")

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

    pheno = ix.phenotypes()[pidx]

    unreachable: list[str] = []

    async def forest_for(idx: int) -> dict | None:
        """None means "no results for this gene"; an outage is recorded apart.

        Collapsing the two made an R2 outage report as "this gene has no BRaVa
        results", which is a claim about biology drawn from a network failure.
        """
        info = ix.gene_info(idx)
        try:
            payload = await client.gene_payload(info["ensg"])
        except client.NotFound:
            return None
        except client.Unavailable:
            unreachable.append(info["gene"])
            return None
        return q.forest(
            payload, pheno, pheno_idx=pidx, mask_idx=mask_idx,
            maf_idx=maf_idx, test=test_name,
        )

    # ---- list form: a compact verdict per gene ----------------------------
    if len(wanted) > 1:
        # Concurrent, but the client coalesces and caps concurrency, so this is
        # still at most one request per distinct gene.
        results = await asyncio.gather(*(forest_for(i) for i in resolved.values()))
        pairs = [
            (ix.gene_info(idx)["gene"], res)
            for (idx, res) in zip(resolved.values(), results)
            if res and res["strata"]
        ]
        missing = [
            g for g, res in zip(resolved, results) if not res or not res["strata"]
        ]
        if not pairs:
            if unreachable:
                return _err(
                    f"BRaVa data is temporarily unreachable for "
                    f"{', '.join(unreachable)}. A fetch failure, not an absence of "
                    "results; retry."
                )
            return _err(
                f"None of those genes has a {MASK_LABEL[mask_idx]} / "
                f"{MAF_LABEL[maf_idx]} result for {pheno['name']}."
            )
        out = {
            "trait": pheno["name"],
            "type": pheno["type"],
            "mask": MASK_LABEL[mask_idx],
            "maf": MAF_LABEL[maf_idx],
            "test": test_name,
            "screened": len(pairs),
            "results": q.replication_summary(pairs),
            "replication_basis": "'concordant' counts the 5 superpopulations whose "
            "effect matches the meta at nominal p<0.05; the verdict is derived by "
            "this server from upstream's numbers, not a published statistic. Call "
            "this tool with a single gene for its full per-ancestry forest.",
            "note": f"{BETA_NOTE} {DISCLAIMER}",
        }
        missing = [g for g in missing if g not in unreachable]
        if missing:
            out["no_result"] = ", ".join(missing)
        if unreachable:
            out["unreachable"] = (
                f"{', '.join(unreachable)}: could not be fetched. Not an absence of "
                "results; retry for these."
            )
        if (mn := _mask_note(mask_idx)):
            out["warning"] = mn
        return _emit(out, narrow="a shorter gene list")

    # ---- single gene: the full forest -------------------------------------
    gidx = next(iter(resolved.values()))
    info = ix.gene_info(gidx)
    result = await forest_for(gidx)
    if result is None:
        if unreachable:
            return _err(
                f"BRaVa data is temporarily unreachable for {info['gene']}. A fetch "
                "failure, not an absence of results; retry."
            )
        return _err(f"{info['gene']} ({info['ensg']}) has no BRaVa results.")
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
        "note": f"{BETA_NOTE} {DISCLAIMER}",
    }
    notes = [n for n in (_ambiguity_note(gene, gidx), _mask_note(mask_idx)) if n]
    if notes:
        out["warning"] = " ".join(notes)
    return toons.dumps(out)


@mcp.tool(annotations=READ_ONLY)
async def top_associations(
    ancestry: str = DEFAULT_ANCESTRY,
    category: str | None = None,
    trait_type: str | None = None,
    genes: str | None = None,
    absent_in: str | None = None,
    mask: str | None = None,
    test: str = DEFAULT_TEST,
    max_p: float = SIG_GENE_CAUCHY,
    limit: int = 25,
    offset: int = 0,
    collapse: bool = True,
    group_by: str | None = None,
) -> str:
    """Strongest rare-variant associations ACROSS traits and genes.

    The cross-cutting view: "the strongest rare-variant signals in cardiovascular
    disease", "which traits does this gene hit at exome-wide significance", "the
    most pleiotropic genes", "screen these candidates", "what is specific to one
    ancestry". Served entirely from the bundled index, with no network access.

    Args:
        ancestry: All (default), EUR, AFR, AMR, EAS, SAS, non_EUR.
        category: Trait category filter, e.g. "Cardiovascular", "Lipids",
            "Endocrine/Metabolic" (see catalog(kind='phenotypes')).
        trait_type: "binary" or "quantitative".
        genes: Restrict to one gene or a comma-separated set ("PCSK9,LDLR,APOB").
            The way to screen a candidate list in a single call instead of one
            call per gene.
        absent_in: Return only findings that clear the threshold in `ancestry`
            but NOT in this stratum, e.g. ancestry="AFR", absent_in="EUR" for
            signals a European-only study would have missed. Ancestry-specific
            effects are the reason the consortium exists, and this is the query
            for them.
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

    gene_idxs = None
    if genes:
        wanted = [g.strip() for g in genes.split(",") if g.strip()]
        resolved = {g: ix.resolve_gene(g) for g in wanted}
        unknown = [g for g, i in resolved.items() if i is None]
        if unknown:
            return _err(
                f"Unknown gene(s): {', '.join(unknown)}. Call search() on each, or "
                "pass Ensembl ids."
            )
        gene_idxs = set(resolved.values())

    absent_idx = None
    if absent_in:
        try:
            absent_idx = ix.resolve_ancestry(absent_in)
        except ValueError as exc:
            return _err(str(exc))
        if absent_idx == anc_idx:
            return _err(
                f"absent_in must differ from ancestry (both are '{absent_in}'). "
                "Use e.g. ancestry='AFR', absent_in='EUR'."
            )

    # Aggregation sorts on p, so it needs the numeric value: formatting first
    # would make it compare "1.00e-08" against "9.99e-300" as strings.
    rows = q.all_results_rows(
        ix.all_results(ANCESTRIES[anc_idx]),
        ix.genes(),
        ix.phenotypes(),
        gene_idxs=gene_idxs,
        mask_idx=mask_idx,
        test_idx=TEST_INDEX[test_name],
        category=category,
        trait_type=trait_type,
        max_p=max_p,
        format_p=group_by is None,
    )
    if group_by is not None and group_by not in ("gene", "trait"):
        return _err(
            f"Unknown group_by '{group_by}'. Valid: 'gene' (rank genes by how many "
            "traits they hit) or 'trait' (rank traits by how many genes reach the "
            "threshold). Omit it for the flat association list."
        )

    contrast_note = None
    if absent_idx is not None:
        other = ANCESTRIES[absent_idx]
        seen = q.significant_pairs(
            ix.all_results(other), max_p, TEST_INDEX[test_name]
        )
        phenos = ix.phenotypes()
        untested = set()
        kept = []
        for row in rows:
            pi = next(
                i for i, ph in enumerate(phenos) if ph["id"] == row["trait_id"]
            )
            gi = ix.resolve_gene(row["ensg"])
            if (gi, pi) in seen:
                continue
            if other not in phenos[pi]["ancestries"]:
                untested.add(row["trait_id"])
            kept.append(row)
        rows = kept
        contrast_note = (
            f"Clears p<{max_p} in {ANCESTRIES[anc_idx]} and not in {other}."
        )
        if untested:
            contrast_note += (
                f" Careful: {other} was never analysed for {', '.join(sorted(untested))}, "
                "so for those the absence is missing data, not a null result."
            )

    total = len(rows)
    # Aggregation needs one row per pair, so the dedupe is not optional there.
    if collapse or group_by:
        rows = q.collapse_best(rows, ("ensg", "trait_id"))
    pairs = len(rows)
    if group_by:
        rows = q.format_pvalues(q.aggregate(rows, group_by))

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
        "contrast": contrast_note,
        "results": rows[offset : offset + limit],
        "note": f"{BETA_NOTE} {DISCLAIMER}",
    }
    if offset + limit < len(rows):
        out["next_offset"] = offset + limit
    return _emit(out, narrow="category=, trait_type=, mask= or a smaller max_p", offset=offset)


@mcp.tool(annotations=READ_ONLY)
async def variants(
    phenotype: str,
    gene: str | None = None,
    ancestry: str = DEFAULT_ANCESTRY,
    chrom: str | None = None,
    max_p: float | None = None,
    limit: int = 25,
    offset: int = 0,
) -> str:
    """Single-variant results for a trait, genome-wide or inside one gene.

    Drops from gene-level burden testing to the individual variants carrying a
    signal. Without `gene` this ranks the whole genome for the trait, which is
    the only way in that does not require knowing the gene first; with `gene` it
    restricts to that gene and additionally reports the per-biobank
    effect-direction tally, the cross-biobank replication evidence.

    Each row links to gnomAD, where population allele frequencies live.

    Args:
        phenotype: Trait id or name.
        gene: Restrict to one gene. Omit for the genome-wide scan.
        ancestry: All (cross-ancestry meta, default) or a specific stratum.
            Only meaningful together with `gene`.
        chrom: Restrict the genome-wide scan to one chromosome ("2", "X").
        max_p: p-value ceiling. The variant-level threshold is 1.82e-8.
        limit: Max rows (default 25).
        offset: Skip this many rows, to page through a long result set.

    Returns: variant (chr-pos-ref-alt), gene, p-value, beta, effect direction,
    the ancestries it was seen in, and a gnomAD link. Within a gene, also the
    95% CI, effective sample size, I-squared, heterogeneity p and the per-biobank
    direction string.
    """
    limit = max(1, min(limit, 200))
    pidx = ix.resolve_phenotype(phenotype)
    if pidx is None:
        return _err(f"Unknown trait '{phenotype}'. Call catalog(kind='phenotypes').")
    pheno = ix.phenotypes()[pidx]

    # ---- genome-wide scan -------------------------------------------------
    if gene is None:
        if ancestry and ancestry.strip().lower() != DEFAULT_ANCESTRY.lower():
            return _err(
                "The genome-wide variant scan is published for the cross-ancestry "
                f"meta only. Pass a gene to see {ancestry} strata, or drop ancestry."
            )
        try:
            payload = await client.variant_overview_payload(pheno["id"])
        except client.NotFound:
            return _err(
                f"No genome-wide variant scan published for {pheno['name']}. "
                f"Pass a gene to get its variants, or use "
                f"phenotype_associations('{pheno['id']}') for gene-level results."
            )
        except client.Unavailable as exc:
            return _err(f"Variant-level data is temporarily unreachable: {exc}")

        rows = vq.overview_rows(
            payload, ix.genes(), pheno["type"], max_p=max_p, chrom=chrom
        )
        out = {
            "trait": pheno["name"],
            "trait_id": pheno["id"],
            "type": pheno["type"],
            "scope": f"genome-wide{f' (chr{chrom})' if chrom else ''}",
            "variant_significance_threshold": SIG_VARIANT,
            "total_matching": len(rows),
            "results": rows[offset : offset + limit],
            "note": f"Upstream thins the null band of this scan, so it ranks real "
            f"signal rather than listing every variant tested. {BETA_NOTE} {DISCLAIMER}",
        }
        if offset + limit < len(rows):
            out["next_offset"] = offset + limit
        return _emit(out, narrow="max_p= or chrom=", offset=offset)

    # ---- within one gene --------------------------------------------------
    gidx = ix.resolve_gene(gene)
    if gidx is None:
        return _err(f"Unknown gene '{gene}'. Call search('{gene}').")

    try:
        anc_idx = ix.resolve_ancestry(ancestry)
    except ValueError as exc:
        return _err(str(exc))
    if anc_idx is None:
        return _err(f"An ancestry is required. Valid: {', '.join(ANCESTRIES)}")

    info = ix.gene_info(gidx)
    split = info["ensg"] in ix.variant_split()
    anc_name = ANCESTRIES[anc_idx]

    try:
        if anc_name == "All":
            payload = await client.gene_variants_payload(
                info["ensg"], pidx if split else None, split
            )
            rows = vq.variant_rows(
                payload, pidx, pheno["type"], max_p=max_p, limit=offset + limit
            )
        else:
            payload = await client.gene_variants_anc_payload(
                info["ensg"], pidx if split else None, split
            )
            rows = vq.ancestry_rows(
                payload,
                pidx,
                anc_name,
                pheno["type"],
                max_p=max_p,
                limit=offset + limit,
                chrom=info["chr"],
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

    out = {
        "gene": info["gene"],
        "ensg": info["ensg"],
        "trait": pheno["name"],
        "type": pheno["type"],
        "ancestry": ANCESTRY_LABEL[anc_name],
        "variant_significance_threshold": SIG_VARIANT,
        "total_matching": len(rows),
        "results": rows[offset : offset + limit],
        "note": f"'biobanks' counts concordant effect directions across contributing "
        f"cohorts; '?' marks a cohort where the variant was absent. "
        f"{BETA_NOTE} {DISCLAIMER}",
    }
    if offset + limit < len(rows):
        out["next_offset"] = offset + limit
    if (amb := _ambiguity_note(gene, gidx)):
        out["warning"] = amb
    return _emit(out, narrow="a smaller max_p", offset=offset)


@mcp.tool(annotations=READ_ONLY)
async def catalog(kind: str = "phenotypes", trait: str | None = None) -> str:
    """What BRaVa covers: traits, contributing biobanks, or the analysis vocabulary.

    Call this before guessing a trait id, or to ground a statement about study
    design (sample sizes, ancestry composition, significance thresholds).

    With kind="biobanks" and a trait, returns which cohorts actually contributed
    to THAT trait's analysis and how many participants and cases each brought.
    That is the question behind any claim of cross-biobank replication: a result
    resting on one large cohort is a different result from one seen in eight.

    Args:
        kind: "phenotypes" (the 44 traits with sample sizes per ancestry),
            "biobanks" (the contributing cohorts), or "vocabulary" (masks, MAF
            cutoffs, tests and significance thresholds).
        trait: With kind="biobanks", restrict to one trait and report each
            cohort's contribution to it, broken down by ancestry.

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
        if trait:
            tidx = ix.resolve_phenotype(trait)
            if tidx is None:
                return _err(
                    f"Unknown trait '{trait}'. Call catalog(kind='phenotypes') for the "
                    "44 traits BRaVa covers."
                )
            pheno = ix.phenotypes()[tidx]
            rows = q.biobank_contributions(ix.pheno_sizes(), ix.biobanks(), pheno["id"])
            if not rows:
                return _err(
                    f"No per-biobank breakdown published for {pheno['name']}. "
                    "catalog(kind='biobanks') gives each cohort's overall size."
                )
            return _emit(
                {
                    "trait": pheno["name"],
                    "trait_id": pheno["id"],
                    "type": pheno["type"],
                    "contributing_biobanks": len(rows),
                    "biobanks": rows,
                    "note": "Per-biobank sizes cover the five superpopulations only; "
                    "the cross-ancestry meta totals in catalog(kind='phenotypes') "
                    "are the authoritative overall N.",
                },
                "biobanks",
            )

        rows = [
            {
                "id": b["id"],
                "name": b["name"],
                "country": b["country"],
                "n": b["sample_size"],
                "sequencing": b["sequencing"],
                "ascertainment": b["ascertainment"],
                "ancestries": ",".join(b["ancestries"]),
                # Which populations a cohort actually holds, and in what numbers.
                # Without it "does any cohort have African-ancestry samples" is
                # unanswerable from a list that only names the strata analysed.
                "ancestry_n": "; ".join(
                    f"{pop} {n:,}"
                    for pop, n in sorted(
                        (b.get("ancestry_n") or {}).items(), key=lambda kv: -kv[1]
                    )
                ),
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
                "data_release": ix.bundle_stamp(),
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
