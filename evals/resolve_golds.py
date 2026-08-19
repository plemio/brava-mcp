"""Recompute the eval gold answers from the RAW upstream files.

Imports nothing from `brava`: it re-reads the published
JSON and re-implements the decoding in a few lines here. If it used the server's
decoder the benchmark would agree with a decoding bug instead of catching it,
and the whole exercise would be circular.

Run it to re-derive the golds after an upstream data freeze, and diff against
evals/questions.json:

    uv run python evals/resolve_golds.py
"""

import collections
import json
import urllib.request
from pathlib import Path

META = "https://nikbaya.github.io/brava_browser/data/meta"
DATA = "https://pub-70f6a636186f47b2a7dbb9547de34be8.r2.dev"
LP_FLOOR = 323.3062153431158


def get(url: str):
    req = urllib.request.Request(url, headers={"User-Agent": "brava-mcp-eval-golds"})
    return json.load(urllib.request.urlopen(req, timeout=120))


def p_of(lp):
    return None if lp is None else (0.0 if lp >= LP_FLOOR else 10.0**-lp)


def main() -> None:
    genes, phen = get(f"{META}/genes.json"), get(f"{META}/phenotypes.json")["phenotypes"]
    allr = get(f"{META}/all_results.All.json")
    sym = lambda gi: genes["symbols"][gi] or genes["ids"][gi]
    pidx = lambda pid: next(i for i, p in enumerate(phen) if p["id"] == pid)
    LDLC, T2D = pidx("LDLC"), pidx("T2Diab")

    # SKAT-O rows, best p per (gene, trait) pair.
    best: dict[tuple[int, int], float] = {}
    for i in range(allr["n"]):
        if allr["test_idx"][i] != 2:
            continue
        key = (allr["gene_idx"][i], allr["pheno_idx"][i])
        p = p_of(allr["lp"][i])
        if key not in best or p < best[key]:
            best[key] = p
    ranked = sorted(best.items(), key=lambda kv: (kv[1], sym(kv[0][0])))

    gene = get(f"{DATA}/gene/ENSG00000169174.json")
    cell = lambda ph, mk, mf, an: next(
        i for i in range(gene["n"])
        if gene["pheno"][i] == ph and gene["mask"][i] == mk
        and gene["maf"][i] == mf and gene["anc"][i] == an
    )
    meta_row = cell(LDLC, 0, 0, 0)
    sign = 1 if gene["beta"][meta_row] > 0 else -1
    concordant = sum(
        1 for i in range(gene["n"])
        if gene["pheno"][i] == LDLC and gene["mask"][i] == 0 and gene["maf"][i] == 0
        and gene["anc"][i] in (1, 2, 3, 4, 5) and gene["beta"][i] is not None
        and (1 if gene["beta"][i] > 0 else -1) == sign
        and gene["lp_skato"][i] is not None and p_of(gene["lp_skato"][i]) < 0.05
    )

    variants = get(f"{DATA}/v2/variant/gene/ENSG00000169174.json")
    sl = variants["by_pheno"][str(LDLC)]
    k = max(range(len(sl["idx"])), key=lambda j: sl["lp"][j] or -1)
    vi = sl["idx"][k]

    afib = next(p for p in phen if p["id"] == "AFib")["n"]["All"]
    plei = collections.Counter(gi for (gi, _), p in best.items() if p < 1.39e-7)
    syn_p = p_of(gene["lp_skato"][cell(LDLC, 3, 0, 0)])

    # --- multi-entity questions (q11-q14) --------------------------------
    sizes = get(f"{META}/pheno_sizes.json")
    banks = {b["id"]: b["name"] for b in get(f"{META}/biobanks.json")["biobanks"]}
    per_bank: dict[str, int] = {}
    for pop_rows in (sizes.get("T2Diab") or {}).values():
        for entry in pop_rows:
            per_bank[entry["id"]] = per_bank.get(entry["id"], 0) + (entry.get("case") or 0)
    top_bank = max(per_bank.items(), key=lambda kv: kv[1])

    # Candidate screen. EVERY candidate in the question is computed, including
    # the one expected to clear: hardcoding it as a fallback would mean the gold
    # no longer depends on the data, which is the whole point of this file.
    gene_p = {}
    for symbol, ensg in (
        ("PCSK9", "ENSG00000169174"),
        ("ACAN", "ENSG00000157766"),
        ("TTN", "ENSG00000155657"),
    ):
        payload = get(f"{DATA}/gene/{ensg}.json")
        # NOT named `best`: that binds the (gene, trait) -> p map above, and
        # shadowing it here would silently corrupt every other gold.
        strongest = min(
            (p_of(payload["lp_skato"][i]) for i in range(payload["n"])
             if payload["pheno"][i] == LDLC and payload["anc"][i] == 0
             and payload["lp_skato"][i] is not None),
            default=None,
        )
        gene_p[symbol] = strongest
    screened = [s for s, v in gene_p.items() if v is not None and v < 2.5e-6]

    # Replication: same-direction and nominally-significant counts per gene.
    def concordance(ensg: str) -> tuple[int, int]:
        payload = get(f"{DATA}/gene/{ensg}.json")
        rows = [
            i for i in range(payload["n"])
            if payload["pheno"][i] == LDLC and payload["mask"][i] == 4
            and payload["maf"][i] == 0
        ]
        meta_i = next(i for i in rows if payload["anc"][i] == 0)
        sign = 1 if payload["beta"][meta_i] > 0 else -1
        strata = [i for i in rows if payload["anc"][i] in (1, 2, 3, 4, 5)]
        same = [i for i in strata if payload["beta"][i] is not None
                and (1 if payload["beta"][i] > 0 else -1) == sign]
        sig = [i for i in same if payload["lp_skato"][i] is not None
               and p_of(payload["lp_skato"][i]) < 0.05]
        return len(sig), len(strata)

    # All five genes the question names, not a sample of them.
    rep = {s: concordance(e) for s, e in (
        ("PCSK9", "ENSG00000169174"),
        ("LDLR", "ENSG00000130164"),
        ("APOB", "ENSG00000084674"),
        ("ANGPTL3", "ENSG00000132855"),
        ("ABCG5", "ENSG00000138075"),
    )}

    # Ancestry contrast, counted as DISTINCT gene-trait pairs, which is what the
    # question asks and what the default (collapsed) response reports.
    afr = get(f"{META}/all_results.AFR.json")
    eur = get(f"{META}/all_results.EUR.json")
    eur_pairs = {
        (eur["gene_idx"][i], eur["pheno_idx"][i])
        for i in range(eur["n"])
        if eur["test_idx"][i] == 2 and p_of(eur["lp"][i]) < 2.5e-6
    }
    afr_pairs = {
        (afr["gene_idx"][i], afr["pheno_idx"][i])
        for i in range(afr["n"])
        if afr["test_idx"][i] == 2 and p_of(afr["lp"][i]) < 2.5e-6
    }
    contrast = len(afr_pairs - eur_pairs)

    golds = {
        "q01": phen[gene["pheno"][meta_row]]["name"],
        "q02": "lowers" if gene["beta"][meta_row] < 0 else "raises",
        "q03": f"{concordant}/5",
        "q04": sym(min((p, gi) for (gi, pi), p in best.items() if pi == T2D)[1]),
        "q05": str(len({gi for (gi, pi), p in best.items() if pi == LDLC and p < 2.5e-6})),
        "q06": sym(plei.most_common(1)[0][0]),
        "q07": f"{afib['n']}/{afib['case']}",
        "q08": f"{sym(ranked[29][0][0])} / {phen[ranked[29][0][1]]['name']}",
        "q09": f'{variants["chr"]}-{variants["pos"][vi]}-{variants["ref"][vi]}-{variants["alt"][vi]}',
        "q10": "no" if syn_p >= 2.5e-6 else "yes",
        "q11": f"{banks[top_bank[0]]} / {top_bank[1]}",
        "q12": " ".join(sorted(screened)),
        "q13": "; ".join(
            f"{g} {rep[g][0]}/{rep[g][1]}"
            for g in ("PCSK9", "LDLR", "APOB", "ANGPTL3", "ABCG5")
        ),
        "q14": str(contrast),
    }

    committed = {q["id"]: q["answer"] for q in json.loads(
        (Path(__file__).parent / "questions.json").read_text())}
    drift = 0
    for qid, gold in golds.items():
        same = committed.get(qid) == gold
        drift += not same
        print(f"{qid}  {'ok ' if same else 'DRIFT'}  computed={gold!r}"
              + ("" if same else f"  committed={committed.get(qid)!r}"))
    print(f"\n{'no drift' if not drift else f'{drift} gold(s) changed upstream'}")
    raise SystemExit(1 if drift else 0)


if __name__ == "__main__":
    main()
