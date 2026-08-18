# BRaVa MCP

An MCP server for the [Biobank Rare Variant Analysis (BRaVa)](https://brava-genetics.github.io/BRaVa/)
consortium's association results: rare **coding**-variant, **gene-based** tests
meta-analysed across ~1.2M individuals from 10 global biobanks, 44 harmonised
traits, 7 ancestry strata.

It exists because no other MCP server covers rare-variant gene-based association
results. gnomAD gives you allele frequencies and constraint; the GWAS Catalog
gives you published common-variant associations; BRaVa gives you what rare coding
variation in a gene does to a trait, plus something no single-programme resource
can: **whether that finding replicates across ancestries and biobanks**.

Summary statistics only. **Not for clinical use.**

## Tools

| Tool | Question |
|---|---|
| `search` | Resolve a gene symbol / Ensembl id / trait name to identifiers |
| `gene_associations` | Phenome-wide scan: which traits is this gene associated with? |
| `phenotype_associations` | Which genes carry rare-variant signal for this trait? |
| `gene_phenotype_detail` | Does this hit replicate across ancestries and biobanks? |
| `top_associations` | Everything cross-cutting, in one call: strongest signals overall, pleiotropy (`group_by`), candidate-list screening (`genes`), and ancestry-specific findings (`absent_in`) |
| `variants` | Single variants, genome-wide for a trait or inside one gene, with per-biobank concordance |
| `catalog` | The 44 traits, the analysis vocabulary, and which cohorts carried a given trait (`trait`) |

### The queries that are one call rather than sixteen

The per-file layout upstream answers "this gene" and "this trait" well and
anything cross-cutting badly, so those compositions live in the server:

```
top_associations(group_by="gene")                  # most pleiotropic genes
top_associations(genes="PCSK9,LDLR,APOB")          # screen a candidate list
top_associations(ancestry="AFR", absent_in="EUR")  # what a EUR-only study missed
catalog(kind="biobanks", trait="T2Diab")           # who contributed, and how much
variants("LDLC")                                   # strongest variants, no gene needed
```

The ancestry contrast is the one to reach for when the question is why the
consortium exists: 27 gene-trait findings clear the gene-level threshold in AFR
and not in EUR. Where a stratum was never analysed for a trait, the response says
so rather than passing missing data off as a null result.

Every tool is read-only and capped at a 25,000-character response. The five that
return ranked rows page with `offset`/`next_offset`; `search` and `catalog` return
bounded sets and do not. When a result is cut, the note names the parameters that
would narrow it and the offset that continues it, rather than just truncating.

## Data source and traffic discipline

All data comes from the public files behind the
[BRaVa browser](https://nikbaya.github.io/brava_browser/)
([source](https://github.com/nikbaya/brava_browser), MIT). Those files sit on a
personal Cloudflare R2 free tier with hard monthly ceilings, so this server is
built to stay far below them:

* **The metadata indexes are committed here** (`data/meta/`), exactly as the
  browser commits them. Because `all_results.{ANC}.json` already carries every
  row clearing the suggestive threshold across all traits and genes, `search`,
  `catalog`, `top_associations` and the default path of
  `phenotype_associations` issue **no outbound request at all**.
* **Everything fetched is cached on disk forever.** The gene-level data is
  immutable per release, so a given gene is downloaded once and never again.
* `tests/test_traffic.py` asserts both of the above, so the commitment is
  checked rather than merely stated. See
  [nikbaya/brava_browser#1](https://github.com/nikbaya/brava_browser/issues/1)
  for the conversation with the upstream author.

Point `BRAVA_DATA_BASE_URL` / `BRAVA_VARIANT_BASE_URL` at a mirror to move off
the upstream bucket entirely, with no code change.

## Running it

```bash
make sync                       # install
make test                       # offline suite
make test-all                   # + live-data checks (wire contract, known biology)
make eval                       # 10 benchmark questions with fixed gold answers
make serve                      # HTTP daemon on :3163
uv run python server.py         # stdio
```

`evals/questions.json` holds ten questions whose answers were resolved directly
from the raw upstream files, independently of this server, so the benchmark
cannot agree with a decoding bug. `evals/selfcheck.py` walks each one through the
tools: currently 10/10, a median of one tool call per question, and zero outbound
HTTP requests for the whole set. It fixes the tool path rather than letting a
model pick it, so it proves the questions are answerable and at what cost, not
that a model finds its way. The model-in-the-loop half is still missing.

| Variable | Default | Purpose |
|---|---|---|
| `MCP_TRANSPORT` | `stdio` | `http` for the shared daemon |
| `MCP_PORT` | `3163` | daemon port |
| `BRAVA_DATA_BASE_URL` | upstream R2 | gene/phenotype files |
| `BRAVA_VARIANT_BASE_URL` | `$DATA/v2` | variant-level files |
| `BRAVA_CACHE_DIR` | `~/.cache/brava-mcp` | disk cache |
| `BRAVA_CACHE_TTL` | `0` (never) | cache expiry, seconds |

## Reading the results

* **beta > 0** increases risk (binary traits) or the trait value (quantitative);
  **beta < 0** decreases it. `beta`/`SE` always come from the inverse-variance-
  weighted **Burden** meta-analysis, even when the reported p-value is SKAT or
  SKAT-O.
* **SKAT-O** is the primary omnibus test and drives significance calls. **Burden**
  is most powerful when a gene's variants point the same way; **SKAT** when they
  are mixed.
* The **`synonymous` mask is a calibration control**. A significant synonymous
  result indicates residual test inflation, not biology.
* Significance thresholds from the flagship paper: gene × mask Bonferroni
  1.39e-7, gene-level Cauchy 2.5e-6, variant-level 1.82e-8.
* BRaVa carries **no allele frequencies** and **no common-variant GWAS**. Variant
  rows link out to gnomAD for the former.

## Architecture

`brava/constants.py`, `query.py` and `variants.py` are pure stdlib: the wire
contract and all decoding, filtering and ranking live there, so the part most
likely to be wrong is unit-testable without a network. `index.py` reads the
bundled metadata. `client.py` is the only module that touches the network.
`server.py` is thin FastMCP wrappers.

The published payloads are columnar: parallel arrays of integer indices into
canonical lists. Upstream documents that encoding as append-only, and
`tests/test_wire_contract.py` re-reads the live data and fails if an index we
would decode falls outside our lists, so a silent relabelling becomes a test
failure rather than a wrong answer.

## Citation

Palmer, Hill, Hodgson, et al. *Rare variant association analyses across 10 global
biobanks*. medRxiv (2026).
[doi:10.64898/2026.05.21.26353759](https://www.medrxiv.org/content/10.64898/2026.05.21.26353759v1.full)

The bundled metadata in `data/meta/` is derived from that release via the BRaVa
browser's ETL, and is redistributed here under the browser's MIT licence.

## Licence

MIT.
