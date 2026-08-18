"""Walk each eval question through the tools and check the evidence comes out.

This is not the full evaluation, and the distinction matters. It checks each
question's `evidence`, the values the tools must RETURN for the answer to be
derivable, never the `answer` itself: several answers are conclusions ("does it
clear the threshold?" -> "no") that no string match can verify, only a model
can. So this proves the data is reachable and at what call cost; whether a model
reaches the right conclusion from it is the other half, and it needs the
model-in-the-loop runner.

A server that answers correctly in 25 calls is a bad server, so the call count is
the number to watch, not just the hit rate.

The LLM-in-the-loop runner (thread_ephemeral against the local engine) is the
missing half; see the skill's references/evaluation.md.

    uv run python evals/selfcheck.py
"""

import asyncio
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import server
from brava import client

QUESTIONS = json.loads((Path(__file__).parent / "questions.json").read_text())


async def paths() -> dict[str, list]:
    """One entry per question: the calls a competent agent would make.

    Most questions are now a single SQL statement. That is the point of the
    rewrite: the ones that used to need a bespoke tool parameter (pleiotropy,
    candidate screening, ancestry contrast, cohort contribution) are clauses.
    """
    return {
        "q01": [server.query(
            "SELECT trait, p_skato FROM results WHERE ensg='ENSG00000169174' "
            "AND ancestry='All' AND mask='pLoF' AND maf='<0.1%' ORDER BY p_skato LIMIT 1")],
        "q02": [server.query(
            "SELECT beta FROM results WHERE gene='PCSK9' AND trait_id='LDLC' "
            "AND ancestry='All' AND mask='pLoF' AND maf='<0.1%'")],
        "q03": [server.gene_phenotype_detail("PCSK9", "LDLC", mask="pLoF", maf="<0.1%")],
        "q04": [server.query(
            "SELECT gene, min(p_skato) p FROM results WHERE trait_id='T2Diab' "
            "AND ancestry='All' AND mask<>'synonymous' GROUP BY gene ORDER BY p LIMIT 1")],
        "q05": [server.query(
            "SELECT count(DISTINCT gene) n FROM results WHERE trait_id='LDLC' "
            "AND ancestry='All' AND p_skato < 2.5e-6")],
        "q06": [server.query(
            "SELECT gene, count(DISTINCT trait) traits FROM results WHERE ancestry='All' "
            "AND p_skato < 1.39e-7 GROUP BY gene ORDER BY traits DESC LIMIT 1")],
        "q07": [server.query(
            "SELECT n, cases FROM pheno_sizes WHERE trait_id='AFib' AND ancestry='All'")],
        "q08": [server.query(
            "SELECT gene, trait, p FROM (SELECT gene, trait, min(p_skato) p FROM results "
            "WHERE ancestry='All' GROUP BY gene, trait) ORDER BY p, gene LIMIT 1 OFFSET 29")],
        "q09": [server.variants("LDLC", gene="PCSK9", limit=1)],
        "q10": [server.query(
            "SELECT p_skato FROM results WHERE gene='PCSK9' AND trait_id='LDLC' "
            "AND ancestry='All' AND mask='synonymous' AND maf='<0.1%'"), server.schema()],
        "q11": [server.query(
            "SELECT biobank, sum(cases) c FROM biobank_sizes WHERE trait_id='T2Diab' "
            "GROUP BY biobank ORDER BY c DESC LIMIT 1")],
        "q12": [server.query(
            "SELECT gene, min(p_skato) p FROM results WHERE ancestry='All' "
            "AND trait_id='LDLC' AND gene IN ('PCSK9','ACAN','TTN') GROUP BY gene ORDER BY p")],
        "q13": [server.gene_phenotype_detail("PCSK9,LDLR,APOB,ANGPTL3,ABCG5", "LDLC")],
        "q14": [server.query(
            "SELECT count(*) n FROM (SELECT DISTINCT a.gene_idx, a.pheno FROM results a "
            "WHERE a.ancestry='AFR' AND a.p_skato < 2.5e-6 AND NOT EXISTS (SELECT 1 FROM "
            "results e WHERE e.ancestry='EUR' AND e.gene_idx=a.gene_idx AND e.pheno=a.pheno "
            "AND e.p_skato < 2.5e-6))")],
    }


def found(gold: str, blob: str) -> bool:
    """Does the value appear in the response, as a value and not inside a word?

    Word-anchored, because a plain substring test cannot fail: "GCK" would be
    satisfied by "GCKR", "27" by the 27 inside "2.27e-05", and "no" by the "no"
    inside "note:". The boundary is widened to exclude "." as well as word
    characters, so a number cannot match part of a longer decimal.

    Deliberately NOT restricted to data lines: `warning:` carries the calibration
    -control flag, which is itself evidence for one of the questions.
    """
    # A number must not match inside a longer one ("27" in "2.27e-05"), so digits
    # get a boundary that also excludes ".". Text must not carry that exclusion:
    # it would refuse "not in EUR" at the end of a sentence, where the next
    # character is a full stop.
    numeric = re.fullmatch(r"[\d.,<>e+-]+", gold) is not None
    right = r"(?![\w.])" if numeric else r"(?!\w)"
    return re.search(rf"(?<![\w.]){re.escape(gold)}{right}", blob) is not None


async def main() -> int:
    client.reset_counters()
    plans = await paths()
    rows, failures = [], 0

    for q in QUESTIONS:
        before = client.outbound_requests()
        calls = plans[q["id"]]
        blob = "".join(await asyncio.gather(*calls))
        ok = all(found(e, blob) for e in q["evidence"])
        failures += not ok
        rows.append(
            (q["id"], "PASS" if ok else "FAIL", len(calls),
             client.outbound_requests() - before, len(blob),
             q["answer"], " + ".join(q["evidence"]))
        )

    print(f"{'id':5} {'':4} {'calls':>5} {'http':>5} {'chars':>7}  {'answer':28} evidence")
    for r in rows:
        print(f"{r[0]:5} {r[1]:4} {r[2]:>5} {r[3]:>5} {r[4]:>7}  {r[5]:28} {r[6]}")
    print()
    print(f"accuracy {len(rows) - failures}/{len(rows)}"
          f" | median calls {sorted(r[2] for r in rows)[len(rows) // 2]}"
          f" | total http {sum(r[3] for r in rows)}"
          f" | total chars {sum(r[4] for r in rows)}")
    await client.close()
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(asyncio.run(main()))
