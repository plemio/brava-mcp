# BRaVa MCP

An MCP server for the [Biobank Rare Variant Analysis (BRaVa)](https://brava-genetics.github.io/BRaVa/)
consortium's association results: rare **coding**-variant, **gene-based** tests
meta-analysed across ~1.2M individuals from 10 global biobanks, 44 harmonised
traits, 7 ancestry strata.

Summary statistics only. **Not for clinical use.**

## SQL over the whole table

`query` runs read-only SQL against all 61,791,444 gene-level rows: every gene x
trait x variant-mask x MAF-cutoff x ancestry cell, with the Burden, SKAT and
SKAT-O p-values, the effect size and its standard error, and the cross-cohort
heterogeneity test. The database is local, so a query costs no network.

```sql
-- what does this gene do
SELECT trait, mask, p_skato, beta FROM results
WHERE gene='PCSK9' AND ancestry='All' AND mask<>'synonymous'
ORDER BY p_skato LIMIT 20

-- most pleiotropic genes
SELECT gene, count(DISTINCT trait) traits FROM results
WHERE ancestry='All' AND p_skato < 1.39e-7 GROUP BY gene ORDER BY traits DESC

-- what a European-only study would have missed
SELECT a.gene, a.trait, a.p_skato FROM results a
WHERE a.ancestry='AFR' AND a.p_skato < 2.5e-6 AND NOT EXISTS (
  SELECT 1 FROM results e WHERE e.ancestry='EUR'
  AND e.gene_idx=a.gene_idx AND e.pheno=a.pheno AND e.p_skato < 2.5e-6)
```

## Tools

| Tool | What it is for |
|---|---|
| `query` | Read-only SQL over the whole gene-level table |
| `schema` | Tables, columns, runnable recipes, and the traps that make a valid query scientifically wrong. **Read this first** |
| `gene_phenotype_detail` | Cross-ancestry replication for a gene-trait pair, or a screen over a hit list |
| `variants` | Single-variant results, genome-wide for a trait or inside one gene |

Why the other three exist.

`gene_phenotype_detail` computes the concordance count over the five
superpopulations only, excluding `All` and `non_EUR`, which pool the same
individuals. A query that aggregates over every ancestry double-counts.

`variants` fetches over HTTP. The variant-level release is a separate upstream
format of 3.09 GiB across ~176,000 objects, against 1.19 GiB for the gene-level
data, and it is not included in the local database.

`schema()` returns the tables and columns, runnable query templates, and ten
pitfalls: effect sizes belong to a different test than the p-value beside them,
one mask is a calibration control rather than a biological category, ancestry
strata overlap, a p-value of exactly zero is the strongest result rather than a
missing one, and six more. Read it before writing SQL.

## The database

873 MB, published as a
[release asset](https://github.com/plemio/brava-mcp/releases/tag/data-v1) and
downloaded once into `~/.cache/brava-mcp/` at first use. Deliberately not
committed: deployments reset the clone on every spawn, so a gigabyte inside it
would be re-fetched forever. Cloning this repo costs 3.9 MB.

Built by `etl/build_db.py` from the 280 published `phenotype/{P}.{ANC}.json`
files. Those carry the same data as the 19,541 per-gene files, so the pivot is
chosen for politeness: 280 requests against 19,541 for identical coverage, once,
rather than one per gene consulted forever. Class B operations are the scarce
resource on the upstream free tier; egress is free on R2.

Sorted on the low-cardinality key columns and built without ART indexes: 2.49 GB
with indexes, 1.75 GB without, 0.87 GB sorted. No index is missed, because these
are filtered scans and DuckDB's zonemaps already serve them. Every query above
returns in under 70 ms.

## Running it

```bash
make sync                 # install
make db                   # download the published database (873 MB, once)
make test                 # offline suite
make test-all             # + live-data checks
make eval                 # 14 benchmark questions, answers derived independently
make serve                # HTTP daemon on :3163
uv run python server.py   # stdio
```

Rebuild the database from upstream with `uv run python etl/build_db.py`
(~200 s: 120 s of downloads, 76 s of loading, then the sorted compaction).

| Variable | Default | Purpose |
|---|---|---|
| `MCP_TRANSPORT` | `stdio` | `http` for the shared daemon |
| `MCP_PORT` | `3163` | daemon port |
| `BRAVA_DB_URL` | the release asset | where to fetch the database |
| `BRAVA_DB_PATH` | `~/.cache/brava-mcp/brava.duckdb` | local database |
| `BRAVA_VARIANT_BASE_URL` | upstream R2 | variant-level files |

## Reading the results

* **beta > 0** increases risk (binary traits) or the trait value (quantitative).
  `beta` and `se` always come from the inverse-variance-weighted **Burden**
  meta-analysis, including on rows where you read `p_skato`. There is no SKAT-O
  effect size.
* **SKAT-O** is the primary omnibus test. **Burden** is most powerful when a
  gene's variants point the same way; **SKAT** when they are mixed.
* The **`synonymous` mask is a calibration control**. A significant synonymous
  result indicates residual test inflation, not biology.
* Thresholds from the flagship paper: gene × mask Bonferroni 1.39e-7, gene-level
  Cauchy 2.5e-6, variant-level 1.82e-8.
* BRaVa carries **no allele frequencies** and **no common-variant GWAS**. Variant
  rows link to gnomAD for the former.

`schema()` returns all of this, plus five more traps, alongside the columns.

## Evaluation

`evals/questions.json` holds fourteen questions, all fourteen resolved directly
from the raw upstream files by `evals/resolve_golds.py`, which imports nothing
from `brava`, so the benchmark cannot agree with a decoding bug and doubles as an
upstream-drift detector.

`evals/selfcheck.py` walks each question through the tools: currently 14/14, a
median of one call per question, and zero outbound HTTP requests for the whole
set. It checks each question's `evidence` (the values the tools must return) and
never its `answer`, because several answers are conclusions no string match can
verify. So it proves the data is reachable and at what cost, not that a model
reaches the right conclusion; that half needs a model-in-the-loop runner and is
still missing.

## Traffic

Gene-level questions are local, so they cost the upstream project nothing at all.
Only `variants` fetches, and each file is cached permanently. Building the
database costs 280 requests, once. See
[nikbaya/brava_browser#1](https://github.com/nikbaya/brava_browser/issues/1) for
the conversation with the upstream author.

## Citation

Palmer, Hill, Hodgson, et al. *Rare variant association analyses across 10 global
biobanks*. medRxiv (2026).
[doi:10.64898/2026.05.21.26353759](https://www.medrxiv.org/content/10.64898/2026.05.21.26353759v1.full)

The database is derived from that release via the BRaVa browser's published
files, and is redistributed under the browser's MIT licence.

## Licence

MIT.
