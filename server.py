"""BRaVa MCP Server.

Gene- and variant-level rare coding-variant association results from the Biobank
Rare Variant Analysis (BRaVa) consortium: ~1.2M individuals, 10 global biobanks,
44 harmonised traits, 7 ancestry strata.

The gene-level results are a single flat fact table, and a model already writes
SQL at expert level, so the table is shipped rather than wrapped. An earlier
version of this server exposed seven hand-carved tools whose parameters
(group_by, absent_in, genes=, collapse) were each a SQL clause reimplemented
worse, one review cycle at a time, while every question nobody had anticipated
stayed out of reach. `query` removes that ceiling.

What survives as a tool is what SQL cannot know: the semantic traps that make a
syntactically perfect query scientifically wrong (`schema`), a derived summary
whose correct form is easy to get subtly wrong (`gene_phenotype_detail`), and the
variant-level data, which lives in a separate volatile format upstream and is
still fetched over HTTP (`variants`).

Data source: the public files behind https://nikbaya.github.io/brava_browser/
(MIT, https://github.com/nikbaya/brava_browser). Summary statistics only.
"""

from __future__ import annotations

import asyncio
import os

import toons
from fastmcp import FastMCP

from brava import client, db, index as ix, query as q, traps, variants as vq
from brava.constants import (
    ANCESTRIES,
    ANCESTRY_LABEL,
    BROWSER_URL,
    CALIBRATION_MASK_INDEX,
    DEFAULT_ANCESTRY,
    DEFAULT_TEST,
    MAF_LABEL,
    MASK_LABEL,
    PAPER_URL,
    SIG_GENE_CAUCHY,
    SIG_GENE_MASK_BONFERRONI,
    SIG_SUGGEST,
    SIG_VARIANT,
    TESTS,
    TEST_LP_KEY,
)

mcp = FastMCP("BRaVa")

DISCLAIMER = "Summary statistics from the BRaVa flagship paper. Not for clinical use."
BETA_NOTE = (
    "beta/SE: inverse-variance-weighted Burden meta-analysis, whichever test the "
    "p-value comes from."
)
CHARACTER_LIMIT = 25_000
READ_ONLY = {"readOnlyHint": True, "idempotentHint": True, "openWorldHint": True}


def _err(message: str) -> str:
    """Errors are guidance: say what went wrong AND what to do instead."""
    return toons.dumps({"error": message})


def _mask_note(mask_idx: int | None) -> str | None:
    if mask_idx == CALIBRATION_MASK_INDEX:
        return traps.SYNONYMOUS_IS_A_CONTROL
    return None


def _compact(value):
    """Render floats compactly, everywhere, at the serialisation boundary.

    A p-value of 1.17e-205 written as a literal decimal costs 170 characters for
    one cell. That bug shipped four times in this server, each time in a new row
    builder that forgot a formatting pass, and the fourth was `query`, where the
    columns are whatever SQL returned and no builder could have known which ones
    were p-values. Fixing it per builder was the mistake: it belongs here, on the
    single path every response takes out.
    """
    if isinstance(value, float):
        if value != value or value in (float("inf"), float("-inf")):
            return None
        # Four significant figures: upstream rounds -log10(p) to two decimals,
        # so more digits than this are precision the data does not have.
        return f"{value:.4g}"
    if isinstance(value, dict):
        return {k: _compact(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_compact(v) for v in value]
    return value


def _emit(payload: dict, rows_key: str = "results", narrow: str = "") -> str:
    """Serialise, and if the result busts the budget, drop rows and say how to narrow."""
    payload = _compact(payload)
    out = toons.dumps(payload)
    rows = payload.get(rows_key)
    if len(out) <= CHARACTER_LIMIT or not isinstance(rows, list) or len(rows) <= 1:
        return out

    def render(count: int) -> str:
        trial = {**payload, rows_key: rows[:count]}
        trial["truncated"] = (
            f"Response exceeded the {CHARACTER_LIMIT}-character budget: showing "
            f"{count} of {len(rows)} rows."
            + (f" Narrow with {narrow}." if narrow else "")
        )
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


# ---------------------------------------------------------------------------

@mcp.tool(annotations=READ_ONLY)
async def query(sql: str, max_rows: int = 50) -> str:
    """Run read-only SQL over the whole BRaVa gene-level results table.

    61.8 million rows: every gene x trait x variant-mask x MAF-cutoff x ancestry
    cell, with the Burden, SKAT and SKAT-O p-values, the effect size and its
    standard error, and the cross-cohort heterogeneity test. Local, so a query
    costs no network.

    **Call `schema()` first.** It returns the tables, the columns, worked query
    templates, and the semantic traps that make a syntactically valid query
    scientifically wrong here. Several of them invert the answer rather than
    degrade it.

    Query the `results` view rather than the raw tables: it exposes p-values
    instead of -log10, and labels instead of integer codes, so the obvious query
    is also the correct one.

    Args:
        sql: One read-only statement (SELECT / WITH / DESCRIBE / SHOW / EXPLAIN).
            Combine steps with a CTE rather than sending several statements.
        max_rows: Rows returned (default 50, capped 500). The response is also
            capped at 25,000 characters, so select the columns you need.

    Returns: the result rows as a table, plus the row count and whether it was cut.
    """
    max_rows = max(1, min(max_rows, 500))
    try:
        columns, rows, truncated = db.run(sql, max_rows)
    except db.UnsafeQuery as exc:
        return _err(str(exc))
    except db.DatabaseUnavailable as exc:
        return _err(str(exc))
    except Exception as exc:  # noqa: BLE001 - the message IS the guidance
        return _err(
            f"Query failed: {exc}. Call schema() for the exact table and column "
            "names, and query the `results` view for labelled values."
        )

    payload: dict = {
        "columns": columns,
        "row_count": len(rows),
        "results": [dict(zip(columns, r)) for r in rows],
    }
    if truncated:
        payload["truncated"] = (
            f"More than {max_rows} rows matched. Raise max_rows, add a LIMIT, or "
            "aggregate in the query."
        )
    payload["note"] = f"{BETA_NOTE} {DISCLAIMER}"
    return _emit(payload, narrow="a tighter WHERE, fewer columns, or an aggregate")


@mcp.tool(annotations=READ_ONLY)
async def schema() -> str:
    """The tables, the query templates, and the traps. Read this before querying.

    Returns the shipped tables with their columns and row counts, worked queries
    for the questions people actually ask, the analysis vocabulary (masks, MAF
    cutoffs, tests, significance thresholds), and a list of ways a correct-looking
    query gives a wrong answer on this data. That list is not boilerplate: it
    covers effect sizes that belong to a different test than the p-value beside
    them, a mask that is a calibration control rather than a biological category,
    ancestry strata that overlap, and p-values of exactly zero that mean the most
    significant result rather than a missing one.

    Returns: tables, columns, recipes, vocabulary, thresholds and pitfalls.
    """
    try:
        tables = db.table_summary()
    except db.DatabaseUnavailable as exc:
        return _err(str(exc))
    return toons.dumps(
        _compact({
            "query_this": "results",
            "tables": tables,
            "recipes": [{"question": qn, "sql": sq} for qn, sq in traps.RECIPES],
            "vocabulary": {
                "ancestries": [{"id": a, "meaning": ANCESTRY_LABEL[a]} for a in ANCESTRIES],
                "masks": MASK_LABEL,
                "maf_cutoffs": MAF_LABEL,
                "tests": [
                    {"name": "SKAT-O", "column": "p_skato",
                     "role": "primary omnibus test; drives significance calls"},
                    {"name": "Burden", "column": "p_burden",
                     "role": "the only test with a directional effect size"},
                    {"name": "SKAT", "column": "p_skat",
                     "role": "most sensitive when a gene mixes risk and protective variants"},
                ],
                "thresholds": {
                    "gene_mask_bonferroni": q.fmt_p(SIG_GENE_MASK_BONFERRONI),
                    "gene_level_cauchy": q.fmt_p(SIG_GENE_CAUCHY),
                    "suggestive": q.fmt_p(SIG_SUGGEST),
                    "variant_level": q.fmt_p(SIG_VARIANT),
                },
            },
            "pitfalls": traps.ALL,
            "also": "gene_phenotype_detail() for cross-ancestry replication, "
                    "variants() for single-variant results.",
            "paper": PAPER_URL,
            "browser": BROWSER_URL,
            "data_release": ix.bundle_stamp(),
            "note": DISCLAIMER,
        })
    )


@mcp.tool(annotations=READ_ONLY)
async def gene_phenotype_detail(
    gene: str,
    phenotype: str,
    mask: str = "pLoF | damaging missense",
    maf: str = "<0.1%",
    test: str = DEFAULT_TEST,
) -> str:
    """Does a gene-trait association replicate across ancestries and biobanks?

    BRaVa's distinctive view, and a tool rather than a documented query because
    the concordance count has to exclude the two pooled strata ('All' and
    'non_EUR') that contain the same individuals as the ones being counted. The
    obvious SQL double-counts and looks entirely reasonable.

    Pass a comma-separated list to screen a whole hit list at once: one gene at a
    time costs a call each, the list form returns a verdict per gene. Verdicts
    separate "underpowered" from "discordant", which is the distinction that
    matters when a stratum is fifteen times smaller than another.

    Args:
        gene: Gene symbol or Ensembl id, or a comma-separated list of them.
        phenotype: Trait id or name.
        mask: Variant annotation mask (default "pLoF | damaging missense").
        maf: "<0.1%" (default) or "<0.01%".
        test: Burden, SKAT or SKAT-O (default SKAT-O).

    Returns: for one gene, every ancestry stratum with its sample size, p-value,
    beta, 95% CI and direction, plus heterogeneity and a concordance count. For a
    list, one verdict row per gene. The concordance is DERIVED here, not published.
    """
    wanted = [g.strip() for g in (gene or "").split(",") if g.strip()]
    if not wanted:
        return _err("Pass a gene symbol or Ensembl id, or a comma-separated list.")
    if len(wanted) > 25:
        return _err(
            f"{len(wanted)} genes is more than this screens at once (max 25). Narrow "
            "the list first with query()."
        )
    resolved = {g: ix.resolve_gene(g) for g in wanted}
    unknown = [g for g, i in resolved.items() if i is None]
    if unknown:
        return _err(
            f"Unknown gene(s): {', '.join(unknown)}. Look them up with "
            "query(\"SELECT symbol, ensg FROM genes WHERE symbol ILIKE '%NAME%'\")."
        )
    seen: set[int] = set()
    resolved = {g: i for g, i in resolved.items() if not (i in seen or seen.add(i))}

    pidx = ix.resolve_phenotype(phenotype)
    if pidx is None:
        return _err(
            f"Unknown trait '{phenotype}'. List them with "
            'query("SELECT trait_id, trait, category FROM phenotypes").'
        )
    pheno = ix.phenotypes()[pidx]

    try:
        mask_idx = ix.resolve_mask(mask)
        maf_idx = ix.resolve_maf(maf)
        test_name = ix.resolve_test(test, DEFAULT_TEST)
    except ValueError as exc:
        return _err(str(exc))
    if mask_idx is None or maf_idx is None:
        return _err("A specific mask and MAF cutoff are required for this view.")

    lp_col = TEST_LP_KEY[test_name].replace("lp_", "p_")
    sizes = {a: v.get("n") for a, v in (pheno.get("n") or {}).items()}

    def forest_for(idx: int) -> dict | None:
        rows = [
            {
                "ancestry": r[0], "p": r[1], "p_het": r[2],
                "beta": r[3], "se": r[4],
            }
            for r in db.connect().execute(
                f"SELECT ancestry, {lp_col}, p_het, beta, se FROM results "
                "WHERE gene_idx = ? AND pheno = ? AND mask_idx = ? AND maf_idx = ?",
                [idx, pidx, mask_idx, maf_idx],
            ).fetchall()
        ]
        if not rows:
            return None
        return q.forest_from_rows(rows, pheno["type"], sizes)

    try:
        forests = {g: forest_for(i) for g, i in resolved.items()}
    except db.DatabaseUnavailable as exc:
        return _err(str(exc))

    pairs = [
        (ix.gene_info(resolved[g])["gene"], f)
        for g, f in forests.items() if f and f["strata"]
    ]
    missing = [ix.gene_info(resolved[g])["gene"] for g, f in forests.items() if not f]

    if not pairs:
        return _err(
            f"No {MASK_LABEL[mask_idx]} / {MAF_LABEL[maf_idx]} result for "
            f"{pheno['name']}. Check which cells exist with "
            f"query(\"SELECT DISTINCT mask, maf FROM results WHERE gene='{wanted[0]}'\")."
        )

    if len(resolved) > 1:
        out = {
            "trait": pheno["name"], "type": pheno["type"],
            "mask": MASK_LABEL[mask_idx], "maf": MAF_LABEL[maf_idx], "test": test_name,
            "screened": len(pairs),
            "results": q.replication_summary(pairs),
            "replication_basis": "'concordant' counts the 5 superpopulations whose "
            "effect matches the meta at nominal p<0.05; the verdict is derived by "
            "this server from upstream's numbers, not a published statistic. Call "
            "this with a single gene for its full per-ancestry forest.",
            "note": f"{BETA_NOTE} {DISCLAIMER}",
        }
        if missing:
            out["no_result"] = ", ".join(missing)
        if (mn := _mask_note(mask_idx)):
            out["warning"] = mn
        return _emit(out, narrow="a shorter gene list")

    only = next(iter(resolved))
    info = ix.gene_info(resolved[only])
    result = forests[only]
    out = {
        "gene": info["gene"], "ensg": info["ensg"],
        "trait": pheno["name"], "type": pheno["type"],
        "mask": MASK_LABEL[mask_idx], "maf": MAF_LABEL[maf_idx], "test": test_name,
        "strata": result["strata"],
        "replication": result["concordance"],
        "note": f"{BETA_NOTE} {DISCLAIMER}",
    }
    notes = [n for n in (_ambiguity_note(only, resolved[only]), _mask_note(mask_idx)) if n]
    if notes:
        out["warning"] = " ".join(notes)
    return toons.dumps(_compact(out))


def _ambiguity_note(query_text: str, chosen: int) -> str | None:
    """Warn when a gene symbol maps to more than one Ensembl gene."""
    others = ix.ambiguous_alternatives(query_text, chosen)
    if not others:
        return None
    return (
        f"'{query_text}' maps to more than one Ensembl gene. These results are for "
        f"{ix.gene_info(chosen)['ensg']}; also matching: {', '.join(others)}. Pass an "
        "Ensembl id to select a specific one."
    )


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

    Drops below the gene-level burden tests to the individual variants carrying a
    signal. Without `gene` this ranks the whole genome for the trait; with `gene`
    it restricts to that gene and adds the per-biobank effect-direction tally,
    the cross-biobank replication evidence.

    Still fetched over HTTP rather than shipped in the database: the variant-level
    format is a separate, actively changing upstream release, an order of
    magnitude larger than the gene-level table, and rebuilt often enough that a
    local copy would be stale within the week. Each file is cached permanently
    once fetched.

    Each row links to gnomAD, where population allele frequencies live.

    Args:
        phenotype: Trait id or name.
        gene: Restrict to one gene. Omit for the genome-wide scan.
        ancestry: All (cross-ancestry meta, default) or a specific stratum. Only
            meaningful together with `gene`.
        chrom: Restrict the genome-wide scan to one chromosome ("2", "X").
        max_p: p-value ceiling. The variant-level threshold is 1.82e-8.
        limit: Max rows (default 25).
        offset: Skip this many rows, to page through a long result set.

    Returns: variant (chr-pos-ref-alt), gene, p-value, beta, effect direction, the
    ancestries it was seen in, and a gnomAD link. Within a gene, also the 95% CI,
    effective sample size, I-squared, heterogeneity p and the per-biobank
    direction string.
    """
    limit = max(1, min(limit, 200))
    pidx = ix.resolve_phenotype(phenotype)
    if pidx is None:
        return _err(
            f"Unknown trait '{phenotype}'. List them with "
            'query("SELECT trait_id, trait FROM phenotypes").'
        )
    pheno = ix.phenotypes()[pidx]

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
                f"No genome-wide variant scan published for {pheno['name']}. Pass a "
                "gene for its variants, or use query() for gene-level results."
            )
        except client.Unavailable as exc:
            return _err(f"Variant-level data is temporarily unreachable: {exc}")

        rows = vq.overview_rows(payload, ix.genes(), pheno["type"], max_p=max_p, chrom=chrom)
        out = {
            "trait": pheno["name"], "trait_id": pheno["id"], "type": pheno["type"],
            "scope": f"genome-wide{f' (chr{chrom})' if chrom else ''}",
            "variant_significance_threshold": q.fmt_p(SIG_VARIANT),
            "total_matching": len(rows),
            "results": rows[offset : offset + limit],
            "note": f"Upstream thins the null band of this scan, so it ranks real "
            f"signal rather than listing every variant tested. {BETA_NOTE} {DISCLAIMER}",
        }
        if offset + limit < len(rows):
            out["next_offset"] = offset + limit
        return _emit(out, narrow="max_p= or chrom=")

    gidx = ix.resolve_gene(gene)
    if gidx is None:
        return _err(
            f"Unknown gene '{gene}'. Look it up with "
            "query(\"SELECT symbol, ensg FROM genes WHERE symbol ILIKE '%NAME%'\")."
        )
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
            payload = await client.gene_variants_payload(info["ensg"], pidx if split else None, split)
            rows = vq.variant_rows(payload, pidx, pheno["type"], max_p=max_p, limit=offset + limit)
        else:
            payload = await client.gene_variants_anc_payload(info["ensg"], pidx if split else None, split)
            rows = vq.ancestry_rows(
                payload, pidx, anc_name, pheno["type"],
                max_p=max_p, limit=offset + limit, chrom=info["chr"],
            )
    except client.NotFound:
        return _err(
            f"No variant-level data for {info['gene']}. Gene-level results are "
            f"unaffected: query them with query()."
        )
    except client.Unavailable as exc:
        return _err(
            f"Variant-level data is temporarily unreachable: {exc}. This is a fetch "
            "failure, not an absence of results; gene-level results are unaffected."
        )
    except (KeyError, TypeError, IndexError) as exc:
        # The upstream variant format is still evolving; fail legibly.
        return _err(
            f"Variant-level data for {info['gene']} could not be decoded ({exc}). The "
            "upstream format may have changed; gene-level results are unaffected."
        )

    if not rows:
        return toons.dumps(_compact({
            "gene": info["gene"], "trait": pheno["name"],
            "ancestry": ANCESTRY_LABEL[anc_name], "results": [],
            "note": f"No variant-level results for this gene x trait x ancestry. {DISCLAIMER}",
        }))

    out = {
        "gene": info["gene"], "ensg": info["ensg"],
        "trait": pheno["name"], "type": pheno["type"],
        "ancestry": ANCESTRY_LABEL[anc_name],
        "variant_significance_threshold": q.fmt_p(SIG_VARIANT),
        "total_matching": len(rows),
        "results": rows[offset : offset + limit],
        "note": f"'biobanks' counts concordant effect directions across contributing "
        f"cohorts; '?' marks a cohort where the variant was absent. {BETA_NOTE} {DISCLAIMER}",
    }
    if offset + limit < len(rows):
        out["next_offset"] = offset + limit
    if (amb := _ambiguity_note(gene, gidx)):
        out["warning"] = amb
    return _emit(out, narrow="a smaller max_p")


if __name__ == "__main__":
    print("BRaVa MCP Server started", flush=True)
    # Transport switch: stdio is the per-worker child the engine spawns, http the
    # long-lived shared daemon reached via MCP_SERVICES_BASE_URL. Public read-only
    # data, no credentials, so a shared daemon is tenant-safe.
    transport = os.getenv("MCP_TRANSPORT", "stdio").lower()
    if transport in ("http", "streamable-http"):
        mcp.run(
            transport="http",
            host=os.getenv("MCP_HOST", "127.0.0.1"),
            port=int(os.getenv("MCP_PORT", "3163")),
        )
    else:
        mcp.run()
