# BRaVa MCP

An MCP server for the [Biobank Rare Variant Analysis (BRaVa)](https://brava-genetics.github.io/BRaVa/)
consortium's association results: rare **coding**-variant, **gene-based** tests
meta-analysed across ~1.2M individuals from 10 global biobanks, 44 harmonised
traits, 7 ancestry strata.

Summary statistics only. **Not for clinical use.**

## It ships the table, not a wrapper around it

The gene-level results are a single flat fact table, and a model already writes
SQL at expert level, so `query` hands it over: 61,791,444 rows, locally, no
network. Any question is a query, including the ones a fixed set of tools would
never have anticipated.

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

The three non-query tools cover what SQL cannot.

`gene_phenotype_detail`, because the concordance count must exclude `All` and
`non_EUR`, which pool the same individuals as the strata being counted: the
obvious SQL double-counts and looks entirely reasonable.

`variants`, because the variant-level release is a separate upstream format, an
order of magnitude larger and rebuilt often enough that a local copy would be
stale within the week.

`schema()`, because a syntactically perfect query can still be scientifically
wrong here. Effect sizes belong to a different test than the p-value beside them,
one mask is a calibration control rather than a biological category, ancestry
strata overlap, and a p-value of exactly zero is the strongest result rather than
a missing one. Several of those invert an answer instead of degrading it, which
is why they travel with the columns rather than sitting in a README.

## The database

Acquired the way `ciqual-mcp` acquires its dataset: downloaded once into
`~/.cache/brava-mcp/` at first use, never shipped in the git clone (the clone is
reset on every daemon spawn, so a gigabyte inside it would be re-fetched
forever). 873 MB, published as a
[release asset](https://github.com/plemio/brava-mcp/releases/tag/data-v1).

Built by `etl/build_db.py` from the 280 published `phenotype/{P}.{ANC}.json`
files. Those carry the same data as the 19,541 per-gene files, so the pivot is
chosen for politeness: 280 requests against 19,541 for identical coverage, once,
rather than one per gene consulted forever. Class B operations are the scarce
resource on the upstream free tier; egress is free on R2.

Sorted on the low-cardinality key columns and built without ART indexes: 2.49 GB
with indexes, 1.75 GB without, 0.87 GB sorted. No index is missed, because these
are filtered scans and DuckDB's zonemaps already serve them. Every query above
returns in under 70 ms.

## Where the numbers come from, and how they go stale

Worth knowing before citing anything from here, because the chain has three links
and only the first is the consortium's.

The **results** are the BRaVa consortium's SAIGE-GENE+ meta-analysis output,
distributed as gzipped TSVs on `gs://brava-meta-analysis` (Requester Pays).

The **derived JSON** this database is built from comes from the browser's ETL,
which is not a passthrough. It makes four decisions that shape what you read:

* The Burden class carries **two** rows per (gene, mask, MAF): a Stouffer row
  with a beta and no standard error, and an inverse-variance-weighted row with
  the real SE and the heterogeneity test. The ETL keeps IVW. Taking the other
  would give different betas and no confidence intervals.
* Gene symbols and coordinates are **not in the results** at all; they are joined
  from Ensembl 110 (GRCh38).
* Trait names, categories and the binary/quantitative split are parsed from the
  [BRaVa curation](https://github.com/BRaVa-genetics/BRaVa_curation) repo.
* SAIGE prints `Pvalue=0` where its own tail computation underflows. The ETL
  floors those rather than nulling them, which is why `p = 0` here means the
  strongest result rather than a missing one.

So this database inherits an interpretation, not only a set of numbers.

**That gives staleness two sources, not one.** A new consortium data freeze is the
obvious one and is roughly annual. The other is a fix to the browser's ETL, which
changes the published JSON with no consortium release behind it: that has
happened four times so far, twice for accuracy bugs their own test suite caught.
The second is the one worth watching, because nothing announces it.

`data/meta/BUNDLE.json` records the upstream `Last-Modified` the bundle was built
against, and `evals/resolve_golds.py` re-derives all fourteen benchmark answers
from the live files and reports drift. That catches a change on fourteen values;
watching the upstream `pipeline/` commits catches it properly.

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
