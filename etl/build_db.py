"""Build brava.duckdb: the whole gene-level fact table as one queryable file.

Why a database rather than a set of hand-carved tools. The gene-level results are
a single flat fact table (gene x trait x mask x MAF x ancestry x test), and a
model already writes SQL at expert level. Every tool parameter the previous
design grew, group_by / absent_in / genes= / collapse, is a SQL clause
reimplemented worse and one review cycle at a time; the queries nobody
anticipated stayed out of reach entirely. Shipping the table removes the ceiling.

Built from `phenotype/{P}.{ANC}.json`, not from `gene/{ENSG}.json`. The two carry
the SAME data pivoted differently, but the phenotype pivot is 280 objects against
19,541 for identical coverage. Class B operations are the scarce resource on the
upstream free tier (egress is free on R2), so this is the polite pivot by two
orders of magnitude, and it happens once rather than per gene consulted forever.

    uv run python etl/build_db.py            # -> build/brava.duckdb
    uv run python etl/build_db.py --limit 6  # a small one, for development
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import duckdb

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from brava.constants import ANCESTRIES, MAFS, MASK_LABEL, MASKS, TESTS  # noqa: E402

R2 = os.getenv("BRAVA_DATA_BASE_URL", "https://pub-70f6a636186f47b2a7dbb9547de34be8.r2.dev")
META = Path(__file__).resolve().parent.parent / "data" / "meta"
STAGE = Path(os.getenv("BRAVA_ETL_STAGE", "/tmp/brava-etl"))
OUT = Path(__file__).resolve().parent.parent / "build" / "brava.duckdb"
UA = {"User-Agent": "brava-mcp-etl/0.2 (+https://github.com/plemio/brava-mcp)"}

# The view everyone should actually query: p-values rather than -log10, labels
# rather than integer codes. It exists so that the obvious query is also the
# correct one, and so that nobody has to know the wire encoding to ask a question.
RESULTS_VIEW = """
CREATE VIEW results AS
SELECT g.symbol AS gene, g.ensg, g.chr,
       p.trait, p.trait_id, p.category, p.type AS trait_type,
       a.ancestry, m.mask_label AS mask, f.label AS maf,
       pow(10, -assoc.lp_skato) AS p_skato,
       pow(10, -assoc.lp_burden) AS p_burden,
       pow(10, -assoc.lp_skat)  AS p_skat,
       pow(10, -assoc.lp_het)   AS p_het,
       assoc.beta, assoc.se,
       assoc.gene_idx, assoc.pheno, assoc.anc,
       assoc.mask AS mask_idx, assoc.maf AS maf_idx
FROM assoc
JOIN genes g USING (gene_idx)
JOIN phenotypes p USING (pheno)
JOIN ancestries a USING (anc)
JOIN masks m ON m.mask = assoc.mask
JOIN mafs f ON f.maf = assoc.maf
"""


def fetch(pheno: str, ancestry: str) -> Path | None:
    """Download one (trait x ancestry) slice, or return None if not published."""
    dest = STAGE / f"{pheno}.{ancestry}.json"
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    url = f"{R2}/phenotype/{pheno}.{ancestry}.json"
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=300) as resp:
            dest.write_bytes(resp.read())
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            return None
        raise
    return dest


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=0, help="stop after N slices")
    parser.add_argument("--jobs", type=int, default=4, help="concurrent downloads")
    args = parser.parse_args()

    STAGE.mkdir(parents=True, exist_ok=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)

    phenotypes = json.loads((META / "phenotypes.json").read_text())["phenotypes"]
    genes = json.loads((META / "genes.json").read_text())

    slices = [
        (pi, p["id"], p["ancestries"].index(a) if a in p["ancestries"] else None, a)
        for pi, p in enumerate(phenotypes)
        for a in ANCESTRIES
        if a in p["ancestries"]
    ]
    if args.limit:
        slices = slices[: args.limit]
    print(f"{len(slices)} trait x ancestry slices to fetch", flush=True)

    started = time.time()
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        paths = list(pool.map(lambda s: fetch(s[1], s[3]), slices))
    print(f"downloaded in {time.time() - started:.0f}s", flush=True)

    if OUT.exists():
        OUT.unlink()
    con = duckdb.connect(str(OUT))

    # --- fact table ---------------------------------------------------------
    con.execute("""
        CREATE TABLE assoc (
            gene_idx INTEGER, pheno UTINYINT, anc UTINYINT,
            mask UTINYINT, maf UTINYINT,
            lp_burden REAL, lp_skat REAL, lp_skato REAL, lp_het REAL,
            beta REAL, se REAL
        )
    """)
    started = time.time()
    for (pi, pheno_id, _, anc), path in zip(slices, paths):
        if path is None:
            continue
        ai = ANCESTRIES.index(anc)
        con.execute(f"CREATE OR REPLACE TEMP TABLE t AS SELECT * FROM read_json('{path}')")
        con.execute(f"""
            INSERT INTO assoc
            SELECT unnest(gene_idx)::INTEGER, {pi}::UTINYINT, {ai}::UTINYINT,
                   unnest(mask)::UTINYINT, unnest(maf)::UTINYINT,
                   unnest(lp_burden)::REAL, unnest(lp_skat)::REAL,
                   unnest(lp_skato)::REAL, unnest(lp_het)::REAL,
                   unnest(beta)::REAL, unnest(se)::REAL
            FROM t
        """)
    print(f"loaded in {time.time() - started:.0f}s", flush=True)

    # --- dimensions ---------------------------------------------------------
    # Denormalised into views so a reader never has to know the integer codes:
    # the whole point of shipping SQL is that nobody decodes indices by hand.
    con.execute("CREATE TABLE genes (gene_idx INTEGER, ensg VARCHAR, symbol VARCHAR, chr VARCHAR, start INTEGER, \"end\" INTEGER)")
    con.executemany(
        "INSERT INTO genes VALUES (?,?,?,?,?,?)",
        [(i, genes["ids"][i], genes["symbols"][i] or genes["ids"][i],
          genes["chr"][i], genes["start"][i], genes["end"][i])
         for i in range(len(genes["ids"]))],
    )
    con.execute("CREATE TABLE phenotypes (pheno UTINYINT, trait_id VARCHAR, trait VARCHAR, category VARCHAR, type VARCHAR)")
    con.executemany(
        "INSERT INTO phenotypes VALUES (?,?,?,?,?)",
        [(i, p["id"], p["name"], p["category"], p["type"]) for i, p in enumerate(phenotypes)],
    )
    con.execute("CREATE TABLE ancestries (anc UTINYINT, ancestry VARCHAR)")
    con.executemany("INSERT INTO ancestries VALUES (?,?)", list(enumerate(ANCESTRIES)))
    con.execute("CREATE TABLE masks (mask UTINYINT, mask_label VARCHAR, raw VARCHAR)")
    con.executemany("INSERT INTO masks VALUES (?,?,?)",
                    [(i, MASK_LABEL[i], MASKS[i]) for i in range(len(MASKS))])
    con.execute("CREATE TABLE mafs (maf UTINYINT, cutoff DOUBLE, label VARCHAR)")
    con.executemany("INSERT INTO mafs VALUES (?,?,?)",
                    [(i, MAFS[i], ["<0.1%", "<0.01%"][i]) for i in range(len(MAFS))])

    # Sample sizes, per trait x ancestry and per contributing biobank.
    con.execute("CREATE TABLE pheno_sizes (trait_id VARCHAR, ancestry VARCHAR, n BIGINT, cases BIGINT, ctrl BIGINT)")
    con.executemany(
        "INSERT INTO pheno_sizes VALUES (?,?,?,?,?)",
        [(p["id"], a, v.get("n"), v.get("case"), v.get("ctrl"))
         for p in phenotypes for a, v in (p.get("n") or {}).items()],
    )
    sizes = json.loads((META / "pheno_sizes.json").read_text())
    banks = {b["id"]: b for b in json.loads((META / "biobanks.json").read_text())["biobanks"]}
    con.execute("CREATE TABLE biobank_sizes (trait_id VARCHAR, superpop VARCHAR, biobank_id VARCHAR, biobank VARCHAR, country VARCHAR, n BIGINT, cases BIGINT)")
    con.executemany(
        "INSERT INTO biobank_sizes VALUES (?,?,?,?,?,?,?)",
        [(tid, pop, e["id"], banks.get(e["id"], {}).get("name", e["id"]),
          banks.get(e["id"], {}).get("country", ""), e.get("n"), e.get("case"))
         for tid, pops in sizes.items() for pop, entries in pops.items() for e in entries],
    )
    con.execute("CREATE TABLE biobanks (id VARCHAR, name VARCHAR, country VARCHAR, sample_size BIGINT, sequencing VARCHAR, ascertainment VARCHAR, ancestries VARCHAR)")
    con.executemany(
        "INSERT INTO biobanks VALUES (?,?,?,?,?,?,?)",
        [(b["id"], b["name"], b["country"], b["sample_size"], b["sequencing"],
          b["ascertainment"], ",".join(b["ancestries"])) for b in banks.values()],
    )

    # (the results view is created on the compacted file, see RESULTS_VIEW)
    _unused = """
        CREATE VIEW results AS
        SELECT g.symbol AS gene, g.ensg, g.chr,
               p.trait, p.trait_id, p.category, p.type AS trait_type,
               a.ancestry, m.mask_label AS mask, f.label AS maf,
               pow(10, -assoc.lp_skato) AS p_skato,
               pow(10, -assoc.lp_burden) AS p_burden,
               pow(10, -assoc.lp_skat) AS p_skat,
               pow(10, -assoc.lp_het) AS p_het,
               assoc.beta, assoc.se,
               assoc.gene_idx, assoc.pheno, assoc.anc, assoc.mask AS mask_idx, assoc.maf AS maf_idx
        FROM assoc
        JOIN genes g USING (gene_idx)
        JOIN phenotypes p USING (pheno)
        JOIN ancestries a USING (anc)
        JOIN masks m ON m.mask = assoc.mask
        JOIN mafs f ON f.maf = assoc.maf
    """

    n = con.execute("SELECT count(*) FROM assoc").fetchone()[0]
    con.execute("CHECKPOINT")
    con.close()

    # Compact into a fresh file, sorted. This is not housekeeping: DuckDB never
    # reclaims space in place, and sorting on the low-cardinality key columns
    # lets them run-length encode. Measured on the full build: 2.49 GB with ART
    # indexes, 1.75 GB without them, 0.87 GB sorted. No index is created, and
    # none is missed: every query below stays under 70 ms on 62M rows because
    # they are filtered scans, which zonemaps already serve.
    print("compacting...", flush=True)
    started = time.time()
    tmp = OUT.with_suffix(".compact")
    tmp.unlink(missing_ok=True)
    con = duckdb.connect(str(tmp))
    con.execute(f"ATTACH '{OUT}' AS src (READ_ONLY)")
    con.execute("CREATE TABLE assoc AS SELECT * FROM src.assoc "
                "ORDER BY pheno, anc, gene_idx, mask, maf")
    for table in ("genes", "phenotypes", "ancestries", "masks", "mafs",
                  "pheno_sizes", "biobank_sizes", "biobanks"):
        con.execute(f"CREATE TABLE {table} AS SELECT * FROM src.{table}")
    con.execute("DETACH src")
    con.execute(RESULTS_VIEW)
    con.execute("CHECKPOINT")
    con.close()
    tmp.replace(OUT)
    print(f"compacted in {time.time() - started:.0f}s", flush=True)
    print(f"{n:,} rows -> {OUT} ({OUT.stat().st_size / 1e9:.2f} GB)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
