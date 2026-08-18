"""Pure decoding of the variant-level (v2) payloads. Stdlib only.

The v2 wire format is denser than the gene-level one: coordinates are stored
ONCE per file in a shared table (pos/ref/alt), and every per-phenotype or
per-ancestry slice references them by integer index. Decoding therefore means
joining a slice back onto that table — get it wrong and you silently attribute
a statistic to the wrong variant, which is the worst failure mode available
here. Hence a separate, fully testable module.

This is also the least stable of the three data tiers: upstream rebuilt and
re-uploaded all ~170k variant objects twice on 2026-08-17. Callers must degrade
gracefully rather than assume the shape holds.

Per-variant `ed` is a per-biobank effect-direction string (one character per
contributing biobank: '+', '-', or '?' for absent) — the cross-biobank
replication signal that no single-programme resource can provide.
"""

from __future__ import annotations

from typing import Any

from .constants import ANCESTRIES, NON_EUR_BIT, SUPERPOPS, p_from_lp, tier
from .query import _sig, ci95, effect_label, format_pvalues

GNOMAD_URL = "https://gnomad.broadinstitute.org/variant/{v}?dataset=gnomad_r4"
CLINVAR_URL = "https://www.ncbi.nlm.nih.gov/clinvar/?term={chrom}[chr]+AND+{pos}[chrpos]"


def decode_anc_mask(mask: int | None) -> str:
    """Presence bitmask -> readable ancestry list.

    Bit i is SUPERPOPS[i]; one bit past them flags "reached the pooled non-EUR
    meta only". A variant with none of the superpop bits but the non-EUR bit set
    is NOT 'no ancestry data' — it is 'non-EUR only', and conflating the two
    was a bug upstream fixed on 2026-08-17.
    """
    if mask is None:
        return ""
    present = [SUPERPOPS[i] for i in range(len(SUPERPOPS)) if mask & (1 << i)]
    if present:
        return ",".join(present)
    if mask & NON_EUR_BIT:
        return "non-EUR (pooled only)"
    return ""


def direction_summary(ed: str | None) -> dict[str, Any] | None:
    """Turn a per-biobank direction string ('++?+-') into a concordance tally."""
    if not ed:
        return None
    plus = ed.count("+")
    minus = ed.count("-")
    tested = plus + minus
    if tested == 0:
        return None
    majority = max(plus, minus)
    return {
        "biobanks": f"{majority}/{tested} concordant",
        "directions": ed,
    }


def variant_rows(
    payload: dict,
    pheno_idx: int,
    trait_type: str,
    *,
    max_p: float | None,
    limit: int,
    chrom: str | None = None,
) -> list[dict[str, Any]]:
    """Decode the all-meta slice for one phenotype inside a gene variant file.

    `payload['by_pheno']` is keyed by the STRINGIFIED phenotype index; a missing
    key means this gene has no variant-level data for that trait, which is a
    legitimate answer, not an error.
    """
    slices = payload.get("by_pheno") or {}
    sl = slices.get(str(pheno_idx))
    if not sl:
        return []

    pos, ref, alt = payload.get("pos", []), payload.get("ref", []), payload.get("alt", [])
    chrom = chrom if chrom is not None else (payload.get("chr") or "")
    out: list[dict[str, Any]] = []

    for k, vi in enumerate(sl.get("idx", [])):
        if vi >= len(pos):
            continue  # defensive: a slice must never outrun its coord table
        p = p_from_lp(_at(sl, "lp", k))
        if max_p is not None and p is not None and p > max_p:
            continue
        beta = _at(sl, "beta", k)
        se = _at(sl, "se", k)
        vid = f"{chrom}-{pos[vi]}-{ref[vi]}-{alt[vi]}"

        row: dict[str, Any] = {
            "variant": vid,
            "chr": chrom,
            "pos": pos[vi],
            "ref": ref[vi],
            "alt": alt[vi],
            "p": p,
            "tier": tier(p),
            "beta": _sig(beta),
            "se": _sig(se),
            "effect": effect_label(beta, trait_type),
        }
        ci = ci95(beta, se)
        if ci:
            row["ci95"] = f"{_sig(ci[0])} to {_sig(ci[1])}"
        row["n_eff"] = _at(sl, "ne", k)
        row["cases"] = _at(sl, "nc", k)
        row["i2"] = _at(sl, "i2", k)
        row["het_p"] = p_from_lp(_at(sl, "cq", k))
        row["ancestries"] = decode_anc_mask(_at(sl, "anc_mask", k))
        dirs = direction_summary(_at(sl, "ed", k))
        if dirs:
            row.update(dirs)
        row["gnomad"] = GNOMAD_URL.format(v=vid)
        out.append(row)

    out.sort(key=lambda r: (r["p"] is None, r["p"]))
    return format_pvalues(out[:limit])


def ancestry_rows(
    payload: dict,
    pheno_idx: int,
    ancestry: str,
    trait_type: str,
    *,
    max_p: float | None,
    limit: int,
    chrom: str = "",
) -> list[dict[str, Any]]:
    """Decode one non-meta ancestry stratum from a `.anc.json` payload.

    `by_anc` is keyed by ancestry index then phenotype index, both stringified.
    """
    anc_idx = ANCESTRIES.index(ancestry)
    by_anc = payload.get("by_anc") or {}
    sl = (by_anc.get(str(anc_idx)) or {}).get(str(pheno_idx))
    if not sl:
        return []

    pos, ref, alt = payload.get("pos", []), payload.get("ref", []), payload.get("alt", [])
    out: list[dict[str, Any]] = []

    for k, vi in enumerate(sl.get("idx", [])):
        if vi >= len(pos):
            continue
        p = p_from_lp(_at(sl, "lp", k))
        if max_p is not None and p is not None and p > max_p:
            continue
        beta = _at(sl, "beta", k)
        se = _at(sl, "se", k)
        row: dict[str, Any] = {
            "variant": f"{chrom}-{pos[vi]}-{ref[vi]}-{alt[vi]}",
            "ancestry": ancestry,
            "p": p,
            "tier": tier(p),
            "beta": _sig(beta),
            "se": _sig(se),
            "effect": effect_label(beta, trait_type),
            "n_eff": _at(sl, "ne", k),
            "cases": _at(sl, "nc", k),
            "i2": _at(sl, "i2", k),
            "het_p": p_from_lp(_at(sl, "cq", k)),
        }
        ci = ci95(beta, se)
        if ci:
            row["ci95"] = f"{_sig(ci[0])} to {_sig(ci[1])}"
        out.append(row)

    out.sort(key=lambda r: (r["p"] is None, r["p"]))
    return format_pvalues(out[:limit])


def _at(slice_: dict, key: str, i: int) -> Any:
    """Positional read that tolerates a column being absent or short.

    The v2 format gained columns twice in a single day upstream; a missing
    column must degrade to None, never raise.
    """
    col = slice_.get(key)
    if not isinstance(col, list) or i >= len(col):
        return None
    return col[i]
