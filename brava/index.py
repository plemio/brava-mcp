"""Bundled metadata indexes + entity resolution.

The `meta/*` files are COMMITTED in this repo (data/meta/), exactly as the BRaVa
browser commits them under app/public/data/meta/. That is deliberate, and it is
what keeps this server's outbound traffic near zero: because
`all_results.{ANC}.json` already carries every row clearing the suggestive
cutoff across all traits and genes, search, the catalogue, top-hits and
cross-trait queries are answered entirely from disk, without touching the
upstream bucket at all.

`make refresh-meta` re-downloads them; that is the only maintenance step, and it
is only needed when the consortium ships a new data freeze.
"""

from __future__ import annotations

import json
import os
import re
from functools import lru_cache
from pathlib import Path
from typing import Any

from .constants import (
    ANCESTRIES,
    MAF_LABEL,
    MAFS,
    MASK_LABEL,
    MASKS,
    TESTS,
)

META_DIR = Path(os.getenv("BRAVA_META_DIR", Path(__file__).resolve().parent.parent / "data" / "meta"))

ENSG_RE = re.compile(r"^ENSG\d{11}$", re.IGNORECASE)


class BravaDataError(RuntimeError):
    """Raised when the bundled metadata is missing or unreadable."""


@lru_cache(maxsize=None)
def _load(name: str) -> Any:
    path = META_DIR / name
    try:
        with path.open(encoding="utf-8") as fh:
            return json.load(fh)
    except FileNotFoundError as exc:
        raise BravaDataError(
            f"Bundled metadata {name} is missing from {META_DIR}. "
            "Run `make refresh-meta` to download it."
        ) from exc


def genes() -> dict:
    """Canonical gene table; array position IS the gene_idx used everywhere."""
    return _load("genes.json")


def phenotypes() -> list[dict]:
    return _load("phenotypes.json")["phenotypes"]


def biobanks() -> list[dict]:
    return _load("biobanks.json")["biobanks"]


def pheno_sizes() -> dict:
    return _load("pheno_sizes.json")


def variant_split() -> set[str]:
    """ENSG ids whose variant data is served per-phenotype rather than in one file."""
    return set(_load("variant_split.json")["split"])


def all_results(ancestry: str) -> dict:
    return _load(f"all_results.{ancestry}.json")


# ---------------------------------------------------------------------------
# Gene resolution
# ---------------------------------------------------------------------------

@lru_cache(maxsize=1)
def _gene_lookup() -> tuple[dict[str, int], dict[str, list[int]]]:
    """(ENSG -> idx, UPPERCASE SYMBOL -> [idx]) built once.

    Symbols map to a LIST: a handful of symbols are reused across Ensembl ids,
    and silently keeping the first would answer a different gene than asked.
    """
    g = genes()
    by_ensg = {ensg.upper(): i for i, ensg in enumerate(g["ids"])}
    by_symbol: dict[str, list[int]] = {}
    for i, sym in enumerate(g["symbols"]):
        if sym:
            by_symbol.setdefault(sym.upper(), []).append(i)
    return by_ensg, by_symbol


def resolve_gene(query: str) -> int | None:
    """Exact resolution of an Ensembl id or gene symbol to a gene_idx."""
    q = (query or "").strip().upper()
    if not q:
        return None
    by_ensg, by_symbol = _gene_lookup()
    if q in by_ensg:
        return by_ensg[q]
    hits = by_symbol.get(q)
    return hits[0] if hits else None


def gene_info(idx: int) -> dict[str, Any]:
    g = genes()
    return {
        "gene": g["symbols"][idx] or g["ids"][idx],
        "ensg": g["ids"][idx],
        "chr": g["chr"][idx],
        "start": g["start"][idx],
        "end": g["end"][idx],
    }


def search_genes(query: str, limit: int) -> list[dict[str, Any]]:
    """Rank gene matches: exact ENSG, exact symbol, prefix, then substring.

    Ranking matters more than it looks: without it a substring hit can outrank
    the exact symbol purely because it sits earlier in the array. Upstream hit
    and fixed exactly this ("Rank exact/prefix gene search matches ahead of
    array-order coincidence").
    """
    q = (query or "").strip().upper()
    if not q:
        return []
    g = genes()
    by_ensg, by_symbol = _gene_lookup()
    seen: set[int] = set()
    out: list[dict[str, Any]] = []

    def push(idx: int, why: str) -> None:
        if idx in seen or len(out) >= limit:
            return
        seen.add(idx)
        out.append({**gene_info(idx), "match": why})

    if q in by_ensg:
        push(by_ensg[q], "exact id")
    for idx in by_symbol.get(q, []):
        push(idx, "exact symbol")
    for sym, idxs in by_symbol.items():
        if len(out) >= limit:
            break
        if sym.startswith(q) and sym != q:
            for idx in idxs:
                push(idx, "prefix")
    for sym, idxs in by_symbol.items():
        if len(out) >= limit:
            break
        if q in sym and not sym.startswith(q):
            for idx in idxs:
                push(idx, "substring")
    return out[:limit]


# ---------------------------------------------------------------------------
# Phenotype resolution
# ---------------------------------------------------------------------------

def resolve_phenotype(query: str) -> int | None:
    """Resolve a trait id ('LDLC') or full name ('LDL cholesterol') to pheno_idx."""
    q = (query or "").strip().lower()
    if not q:
        return None
    phenos = phenotypes()
    for i, p in enumerate(phenos):
        if p["id"].lower() == q:
            return i
    for i, p in enumerate(phenos):
        if p["name"].lower() == q:
            return i
    matches = [i for i, p in enumerate(phenos) if q in p["name"].lower()]
    return matches[0] if len(matches) == 1 else None


def search_phenotypes(query: str, limit: int) -> list[dict[str, Any]]:
    q = (query or "").strip().lower()
    if not q:
        return []
    out: list[dict[str, Any]] = []
    for i, p in enumerate(phenotypes()):
        if q == p["id"].lower():
            why = "exact id"
        elif q == p["name"].lower():
            why = "exact name"
        elif q in p["name"].lower() or q in p["category"].lower():
            why = "substring"
        else:
            continue
        n_all = (p.get("n") or {}).get("All", {})
        out.append(
            {
                "trait_id": p["id"],
                "trait": p["name"],
                "category": p["category"],
                "type": p["type"],
                "n": n_all.get("n"),
                "cases": n_all.get("case"),
                "match": why,
            }
        )
    order = {"exact id": 0, "exact name": 1, "substring": 2}
    out.sort(key=lambda r: order[r["match"]])
    return out[:limit]


# ---------------------------------------------------------------------------
# Vocabulary resolution: accept what a human would type
# ---------------------------------------------------------------------------

def resolve_ancestry(value: str | None) -> int | None:
    """Ancestry name -> index. None/'' means 'no filter'."""
    if value is None or value == "":
        return None
    v = value.strip()
    for i, a in enumerate(ANCESTRIES):
        if a.lower() == v.lower() or a.lower().replace("_", "-") == v.lower():
            return i
    raise ValueError(f"Unknown ancestry '{value}'. Valid: {', '.join(ANCESTRIES)}")


def resolve_mask(value: str | None) -> int | None:
    """Accepts the raw Group string or the short label ('pLoF', 'all variants')."""
    if value is None or value == "":
        return None
    v = value.strip().lower()
    for i, raw in enumerate(MASKS):
        if raw.lower() == v or MASK_LABEL[i].lower() == v:
            return i
    raise ValueError(f"Unknown mask '{value}'. Valid: {', '.join(MASK_LABEL)}")


def resolve_maf(value: str | float | None) -> int | None:
    """Accepts 0.001 / '0.001' / '<0.1%' / '0.1%'."""
    if value is None or value == "":
        return None
    if isinstance(value, (int, float)):
        for i, m in enumerate(MAFS):
            if abs(float(value) - m) < 1e-12:
                return i
        raise ValueError(f"Unknown MAF cutoff {value}. Valid: {MAFS}")
    v = str(value).strip().lower().lstrip("<").strip()
    for i, label in enumerate(MAF_LABEL):
        if label.lower() == str(value).strip().lower() or label.lower().lstrip("<") == v:
            return i
    try:
        return resolve_maf(float(v))
    except ValueError:
        raise ValueError(f"Unknown MAF cutoff '{value}'. Valid: {', '.join(MAF_LABEL)}")


def resolve_test(value: str | None, default: str) -> str:
    if value is None or value == "":
        return default
    v = value.strip().lower().replace("_", "-")
    for t in TESTS:
        if t.lower() == v or t.lower().replace("-", "") == v.replace("-", ""):
            return t
    raise ValueError(f"Unknown test '{value}'. Valid: {', '.join(TESTS)}")
