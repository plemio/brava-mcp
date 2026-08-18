"""The local DuckDB: acquisition, read-only access, and guarded queries.

Acquired the way ciqual-mcp acquires its dataset: downloaded once into the user
cache at first use, never shipped in the git clone. The clone is reset on every
daemon spawn, so a gigabyte inside it would be re-fetched forever; a cache entry
beside it is fetched once per machine.

Read-only is enforced twice, because once is an accident away from being zero:
DuckDB opens the file in read_only mode, AND the statement is checked before it
runs. The connection guarantee is the one that actually holds; the statement
check exists to fail with something a model can act on rather than a driver
error.
"""

from __future__ import annotations

import os
import re
import shutil
import threading
import urllib.request
from pathlib import Path
from typing import Any

import duckdb

DB_URL = os.getenv(
    "BRAVA_DB_URL",
    "https://github.com/plemio/brava-mcp/releases/download/data-v1/brava.duckdb",
)
CACHE_DIR = Path(os.getenv("BRAVA_CACHE_DIR", Path.home() / ".cache" / "brava-mcp"))
DB_PATH = Path(os.getenv("BRAVA_DB_PATH", CACHE_DIR / "brava.duckdb"))
USER_AGENT = "brava-mcp/0.2 (+https://github.com/plemio/brava-mcp)"

# Anything that could write, attach another file, read the filesystem, or shell
# out. DuckDB's read_only connection already refuses the writes; this list is
# about answering with guidance instead of a driver exception, and about the
# functions that read outside the database even on a read-only connection.
_FORBIDDEN = re.compile(
    r"\b(attach|detach|copy|export|import|install|load|"
    r"create|insert|update|delete|drop|alter|truncate|"
    r"read_csv|read_json|read_parquet|read_text|read_blob|glob)\b",
    re.IGNORECASE,
)
_STARTS_OK = re.compile(r"^\s*(select|with|describe|show|explain|pragma\s+table_info)\b", re.IGNORECASE)

_con: duckdb.DuckDBPyConnection | None = None
# Serialises acquisition: the readiness banner prints immediately and the daemon
# starts fetching in the background, so the first real query usually finds the
# file already there instead of waiting out an 873 MB download inside a tool
# timeout. Concurrent callers wait on the same download rather than starting
# their own.
_download_lock = threading.Lock()


class DatabaseUnavailable(RuntimeError):
    """The local database is absent and could not be fetched."""


class UnsafeQuery(ValueError):
    """The statement was refused before it ran."""


def ensure_database() -> Path:
    """Return the local database path, downloading it once if needed."""
    if DB_PATH.exists() and DB_PATH.stat().st_size > 0:
        return DB_PATH
    with _download_lock:
        # Re-checked under the lock: whoever held it may have just finished.
        if DB_PATH.exists() and DB_PATH.stat().st_size > 0:
            return DB_PATH
        return _download()


def prefetch() -> threading.Thread | None:
    """Start acquiring the database in the background, if it is not already here.

    Called at daemon startup. The readiness probe matches the banner printed
    before this returns, so the daemon reports healthy while the file lands, and
    a first query arriving mid-download blocks on the same lock rather than
    starting a second one.
    """
    if DB_PATH.exists() and DB_PATH.stat().st_size > 0:
        return None
    thread = threading.Thread(target=_quiet_prefetch, daemon=True, name="brava-db-fetch")
    thread.start()
    return thread


def _quiet_prefetch() -> None:
    try:
        ensure_database()
    except DatabaseUnavailable:
        # Reported properly when a tool is actually called; failing loudly here
        # would only put a traceback in the daemon log at boot.
        pass


def _download() -> Path:
    DB_PATH.parent.mkdir(parents=True, exist_ok=True)
    tmp = DB_PATH.with_suffix(".part")
    try:
        request = urllib.request.Request(DB_URL, headers={"User-Agent": USER_AGENT})
        with urllib.request.urlopen(request, timeout=1800) as resp, tmp.open("wb") as out:
            shutil.copyfileobj(resp, out, length=1 << 20)
        tmp.replace(DB_PATH)
    except Exception as exc:  # noqa: BLE001 - surfaced to the model as guidance
        tmp.unlink(missing_ok=True)
        raise DatabaseUnavailable(
            f"Could not fetch the BRaVa database from {DB_URL}: {exc}. "
            "Set BRAVA_DB_PATH to a local copy, or BRAVA_DB_URL to a mirror."
        ) from exc
    return DB_PATH


def connect() -> duckdb.DuckDBPyConnection:
    global _con
    if _con is None:
        _con = duckdb.connect(str(ensure_database()), read_only=True)
    return _con


def check_statement(sql: str) -> None:
    """Refuse anything that is not a read, with a message that says what to do."""
    if not sql or not sql.strip():
        raise UnsafeQuery("Empty query. Call schema() to see the tables and columns.")
    if ";" in sql.strip().rstrip(";"):
        raise UnsafeQuery(
            "One statement per call. Split them, or combine them with a CTE (WITH ...)."
        )
    if not _STARTS_OK.match(sql):
        raise UnsafeQuery(
            "Only read statements are accepted: SELECT, WITH, DESCRIBE, SHOW, EXPLAIN. "
            "This database is a published snapshot and is never modified."
        )
    if (hit := _FORBIDDEN.search(sql)) is not None:
        raise UnsafeQuery(
            f"'{hit.group(0)}' is not available here: the database is opened read-only "
            "and file access is disabled. Query the shipped tables instead; schema() "
            "lists them."
        )


def run(sql: str, max_rows: int) -> tuple[list[str], list[tuple[Any, ...]], bool]:
    """Execute a checked read. Returns (columns, rows, truncated).

    Fetches one row beyond the cap so truncation is a fact rather than a guess.
    """
    check_statement(sql)
    cursor = connect().execute(sql)
    rows = cursor.fetchmany(max_rows + 1)
    columns = [d[0] for d in cursor.description] if cursor.description else []
    truncated = len(rows) > max_rows
    return columns, rows[:max_rows], truncated


def table_summary() -> list[dict[str, Any]]:
    """Row counts per shipped table, for the schema tool."""
    con = connect()
    names = [r[0] for r in con.execute(
        "SELECT table_name FROM information_schema.tables "
        "WHERE table_schema='main' ORDER BY table_name"
    ).fetchall()]
    out = []
    for name in names:
        count = con.execute(f'SELECT count(*) FROM "{name}"').fetchone()[0]
        cols = con.execute(f'PRAGMA table_info("{name}")').fetchall()
        out.append({
            "table": name,
            "rows": count,
            "columns": ", ".join(f"{c[1]} {c[2]}" for c in cols),
        })
    return out


def close() -> None:
    global _con
    if _con is not None:
        _con.close()
        _con = None
