"""HTTP access to the bulky per-gene / per-phenotype / per-variant payloads.

Everything answerable from the bundled metadata is answered in `index`; this
module is only reached for a SPECIFIC gene's or trait's detail file.

Traffic discipline — this is a commitment we made to the upstream author, whose
data sits on a personal Cloudflare R2 free tier with hard monthly ceilings:

  * every fetched file is cached on disk FOREVER by default. The v1 gene-level
    data is immutable (every object still carries Last-Modified: 12 Aug 2026),
    so a given gene is downloaded once per deployment and never again.
  * `outbound_requests()` counts what actually left the process, so the promise
    is testable rather than merely stated (see tests/test_traffic.py).
  * a semaphore keeps us from ever looking like a scraper.

Pointing BRAVA_DATA_BASE_URL / BRAVA_META_BASE_URL / BRAVA_VARIANT_BASE_URL at a
mirror moves all of this off the upstream bucket with no code change — the same
escape hatch upstream gives its own frontend via VITE_DATA_BASE_URL.
"""

from __future__ import annotations

import asyncio
import hashlib
import json
import os
import ssl
import time
from pathlib import Path
from typing import Any

import aiohttp
import certifi

# Upstream defaults. Both are public and unauthenticated.
DATA_BASE = os.getenv(
    "BRAVA_DATA_BASE_URL", "https://pub-70f6a636186f47b2a7dbb9547de34be8.r2.dev"
).rstrip("/")
VARIANT_BASE = os.getenv("BRAVA_VARIANT_BASE_URL", f"{DATA_BASE}/v2").rstrip("/")

CACHE_DIR = Path(
    os.getenv("BRAVA_CACHE_DIR", Path.home() / ".cache" / "brava-mcp")
)
# 0 (the default) means "never expire" — correct for immutable release data.
CACHE_TTL = int(os.getenv("BRAVA_CACHE_TTL", "0"))

USER_AGENT = os.getenv(
    "BRAVA_USER_AGENT",
    "brava-mcp/0.1 (+https://github.com/plemio/brava-mcp; MCP server for the BRaVa browser data)",
)

SEM = asyncio.Semaphore(4)
SSL_CTX = ssl.create_default_context(cafile=certifi.where())
TIMEOUT = aiohttp.ClientTimeout(total=90, connect=10, sock_read=60)

_session: aiohttp.ClientSession | None = None
_memory: dict[str, Any] = {}
_outbound = 0


class NotFound(LookupError):
    """The upstream bucket has no object at that path (HTTP 404)."""


class Unavailable(RuntimeError):
    """Upstream is unreachable or returned something unusable."""


def outbound_requests() -> int:
    """Number of HTTP requests this process has actually sent upstream."""
    return _outbound


def reset_counters() -> None:
    _outbound_reset()


def _outbound_reset() -> None:
    global _outbound
    _outbound = 0


def _cache_path(url: str) -> Path:
    digest = hashlib.sha256(url.encode()).hexdigest()[:24]
    return CACHE_DIR / f"{digest}.json"


def _read_disk(url: str) -> Any | None:
    path = _cache_path(url)
    if not path.exists():
        return None
    if CACHE_TTL and (time.time() - path.stat().st_mtime) > CACHE_TTL:
        return None
    try:
        with path.open(encoding="utf-8") as fh:
            return json.load(fh)
    except (OSError, json.JSONDecodeError):
        # A truncated cache entry must never masquerade as upstream data.
        path.unlink(missing_ok=True)
        return None


def _write_disk(url: str, payload: Any) -> None:
    try:
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        tmp = _cache_path(url).with_suffix(".tmp")
        with tmp.open("w", encoding="utf-8") as fh:
            json.dump(payload, fh, separators=(",", ":"))
        tmp.replace(_cache_path(url))
    except OSError:
        # A read-only or full cache dir degrades to memory-only, never fatally.
        pass


def _get_session() -> aiohttp.ClientSession:
    global _session
    if _session is None or _session.closed:
        _session = aiohttp.ClientSession(
            connector=aiohttp.TCPConnector(ssl=SSL_CTX, limit=8),
            timeout=TIMEOUT,
            headers={"User-Agent": USER_AGENT, "Accept-Encoding": "gzip"},
        )
    return _session


async def fetch_json(url: str) -> Any:
    """Fetch and cache one JSON payload. Memory, then disk, then the network."""
    global _outbound

    if url in _memory:
        return _memory[url]

    cached = _read_disk(url)
    if cached is not None:
        _memory[url] = cached
        return cached

    try:
        async with SEM:
            _outbound += 1
            async with _get_session().get(url) as resp:
                if resp.status == 404:
                    raise NotFound(url)
                resp.raise_for_status()
                payload = await resp.json(content_type=None)
    except NotFound:
        raise
    except asyncio.TimeoutError as exc:
        raise Unavailable(f"Timed out fetching {url}") from exc
    except aiohttp.ClientError as exc:
        raise Unavailable(f"Could not fetch {url}: {exc}") from exc

    _memory[url] = payload
    _write_disk(url, payload)
    return payload


# ---------------------------------------------------------------------------
# Typed accessors — the only places that know the upstream path layout
# ---------------------------------------------------------------------------

async def gene_payload(ensg: str) -> dict:
    """gene/{ENSG}.json — every result row for one gene."""
    return await fetch_json(f"{DATA_BASE}/gene/{ensg}.json")


async def phenotype_payload(pheno_id: str, ancestry: str) -> dict:
    """phenotype/{P}.{ANC}.json — every gene for one trait x ancestry (~2.3 MB)."""
    return await fetch_json(f"{DATA_BASE}/phenotype/{pheno_id}.{ancestry}.json")


async def gene_variants_payload(ensg: str, pheno_idx: int | None, split: bool) -> dict:
    """v2/variant/gene/{ENSG}[.{pheno_idx}].json — all-meta variants for a gene.

    Oversized genes are served per-phenotype; `split` comes from the bundled
    variant_split.json manifest.
    """
    name = f"{ensg}.{pheno_idx}" if split else ensg
    return await fetch_json(f"{VARIANT_BASE}/variant/gene/{name}.json")


async def gene_variants_anc_payload(ensg: str, pheno_idx: int | None, split: bool) -> dict:
    """v2/variant/gene/{ENSG}[.{pheno_idx}].anc.json — the non-meta strata."""
    name = f"{ensg}.{pheno_idx}" if split else ensg
    return await fetch_json(f"{VARIANT_BASE}/variant/gene/{name}.anc.json")


async def close() -> None:
    global _session
    if _session and not _session.closed:
        await _session.close()
    _session = None
