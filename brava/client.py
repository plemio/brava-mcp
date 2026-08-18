"""HTTP access to the bulky per-gene / per-phenotype / per-variant payloads.

Everything answerable from the bundled metadata is answered in `index`; this
module is only reached for a SPECIFIC gene's or trait's detail file.

Traffic discipline. This is a commitment we made to the upstream author, whose
data sits on a personal Cloudflare R2 free tier with hard monthly ceilings:

  * every fetched file is cached on disk FOREVER by default. The v1 gene-level
    data is immutable (every object still carries Last-Modified: 12 Aug 2026),
    so a given gene is downloaded once per deployment and never again.
  * concurrent fetches of the same URL are COALESCED onto one request. Agents
    routinely issue tool calls in parallel, and without this two calls landing on
    a cold gene both miss the cache and both go out, which quietly breaks "once,
    ever" the moment anything is concurrent.
  * a 404 is cached too. Roughly 500 of the 20,033 Ensembl genes were never
    tested, so an agent sweeping a gene list hits them repeatedly; without a
    negative cache each miss is an unbounded, repeatable request upstream.
  * `outbound_requests()` counts what actually left the process, so the promise
    is testable rather than merely stated (see tests/test_traffic.py).
  * a semaphore keeps us from ever looking like a scraper.

Pointing BRAVA_DATA_BASE_URL / BRAVA_META_BASE_URL / BRAVA_VARIANT_BASE_URL at a
mirror moves all of this off the upstream bucket with no code change, the same
escape hatch upstream gives its own frontend via VITE_DATA_BASE_URL.
"""

from __future__ import annotations

import asyncio
import hashlib
import json
import os
import ssl
import time
from collections import OrderedDict
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
# 0 (the default) means "never expire", correct for immutable release data.
CACHE_TTL = int(os.getenv("BRAVA_CACHE_TTL", "0"))

USER_AGENT = os.getenv(
    "BRAVA_USER_AGENT",
    "brava-mcp/0.1 (+https://github.com/plemio/brava-mcp; MCP server for the BRaVa browser data)",
)

SEM = asyncio.Semaphore(4)
SSL_CTX = ssl.create_default_context(cafile=certifi.where())
TIMEOUT = aiohttp.ClientTimeout(total=90, connect=10, sock_read=60)

_session: aiohttp.ClientSession | None = None
_session_loop: "asyncio.AbstractEventLoop | None" = None
# Bounded: a phenotype payload is ~2.3 MB and this process is a long-lived
# daemon, so an unbounded dict would grow into the PM2 memory ceiling. Eviction
# only costs a disk read, never a request upstream.
_MEMORY_MAX = 24
_memory: "OrderedDict[str, Any]" = OrderedDict()
# URL -> the task already fetching it, so concurrent callers await one request.
_inflight: dict[str, "asyncio.Task[Any]"] = {}
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


def _remember(url: str, payload: Any) -> None:
    _memory[url] = payload
    _memory.move_to_end(url)
    while len(_memory) > _MEMORY_MAX:
        _memory.popitem(last=False)


def _outbound_reset() -> None:
    global _outbound
    _outbound = 0


def _cache_path(url: str) -> Path:
    digest = hashlib.sha256(url.encode()).hexdigest()[:24]
    return CACHE_DIR / f"{digest}.json"


def _missing_path(url: str) -> Path:
    """Marker for a URL upstream has no object at.

    Kept on disk rather than in memory so it also survives the per-worker stdio
    fallback, where the process is torn down after every call and an in-memory
    set would never be consulted twice.
    """
    digest = hashlib.sha256(url.encode()).hexdigest()[:24]
    return CACHE_DIR / f"{digest}.404"


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
    """The shared session, rebuilt if the running event loop changed.

    An aiohttp session binds to the loop that created it; reusing it from a
    different one raises "Event loop is closed" on the first request. The daemon
    only ever has one loop, but the stdio fallback and the test suite do not, and
    a cache of connections is not worth a hard failure.
    """
    global _session, _session_loop
    loop = asyncio.get_running_loop()
    if _session is None or _session.closed or _session_loop is not loop:
        _session = aiohttp.ClientSession(
            connector=aiohttp.TCPConnector(ssl=SSL_CTX, limit=8),
            timeout=TIMEOUT,
            headers={"User-Agent": USER_AGENT, "Accept-Encoding": "gzip"},
        )
        _session_loop = loop
    return _session


async def fetch_json(url: str) -> Any:
    """Fetch one JSON payload. Memory, then disk, then a coalesced request."""
    if url in _memory:
        _memory.move_to_end(url)
        return _memory[url]

    cached = _read_disk(url)
    if cached is not None:
        _remember(url, cached)
        return cached

    if _missing_path(url).exists():
        raise NotFound(url)

    existing = _inflight.get(url)
    if existing is not None:
        return await asyncio.shield(existing)

    task = asyncio.ensure_future(_fetch_uncached(url))
    _inflight[url] = task
    try:
        return await asyncio.shield(task)
    finally:
        # Only the originator clears the slot, and only once the task is done,
        # so a second caller arriving mid-flight still finds it.
        if task.done():
            _inflight.pop(url, None)


async def _fetch_uncached(url: str) -> Any:
    global _outbound

    try:
        async with SEM:
            _outbound += 1
            async with _get_session().get(url) as resp:
                if resp.status == 404:
                    _mark_missing(url)
                    raise NotFound(url)
                resp.raise_for_status()
                payload = await resp.json(content_type=None)
    except NotFound:
        raise
    except asyncio.TimeoutError as exc:
        raise Unavailable(f"Timed out fetching {url}") from exc
    except aiohttp.ClientError as exc:
        raise Unavailable(f"Could not fetch {url}: {exc}") from exc

    _remember(url, payload)
    _write_disk(url, payload)
    return payload


def _mark_missing(url: str) -> None:
    try:
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        _missing_path(url).touch()
    except OSError:
        pass


# ---------------------------------------------------------------------------
# Typed accessors: the only places that know the upstream path layout
# ---------------------------------------------------------------------------

async def gene_payload(ensg: str) -> dict:
    """gene/{ENSG}.json holds every result row for one gene."""
    return await fetch_json(f"{DATA_BASE}/gene/{ensg}.json")


async def phenotype_payload(pheno_id: str, ancestry: str) -> dict:
    """phenotype/{P}.{ANC}.json holds every gene for one trait x ancestry (~2.3 MB)."""
    return await fetch_json(f"{DATA_BASE}/phenotype/{pheno_id}.{ancestry}.json")


async def gene_variants_payload(ensg: str, pheno_idx: int | None, split: bool) -> dict:
    """v2/variant/gene/{ENSG}[.{pheno_idx}].json holds the all-meta variants for a gene.

    Oversized genes are served per-phenotype; `split` comes from the bundled
    variant_split.json manifest.
    """
    name = f"{ensg}.{pheno_idx}" if split else ensg
    return await fetch_json(f"{VARIANT_BASE}/variant/gene/{name}.json")


async def gene_variants_anc_payload(ensg: str, pheno_idx: int | None, split: bool) -> dict:
    """v2/variant/gene/{ENSG}[.{pheno_idx}].anc.json holds the non-meta strata."""
    name = f"{ensg}.{pheno_idx}" if split else ensg
    return await fetch_json(f"{VARIANT_BASE}/variant/gene/{name}.anc.json")


async def close() -> None:
    global _session, _session_loop
    _inflight.clear()
    if _session and not _session.closed:
        await _session.close()
    _session = None
    _session_loop = None
