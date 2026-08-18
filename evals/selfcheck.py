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
    """One entry per question: the tool calls a competent agent would make."""
    return {
        "q01": [server.gene_associations("ENSG00000169174", mask="pLoF", maf="<0.1%", limit=1)],
        "q02": [server.gene_associations("PCSK9", mask="pLoF", maf="<0.1%", limit=1)],
        "q03": [server.gene_phenotype_detail("PCSK9", "LDLC", mask="pLoF", maf="<0.1%")],
        "q04": [server.phenotype_associations("Type 2 diabetes", limit=1)],
        "q05": [server.phenotype_associations("LDLC", max_p=2.5e-6, limit=1)],
        "q06": [server.top_associations(max_p=1.39e-7, group_by="gene", limit=1)],
        "q07": [server.catalog("phenotypes")],
        "q08": [server.top_associations(max_p=1e-4, limit=1, offset=29)],
        "q09": [server.variants("LDLC", gene="PCSK9", limit=1)],
        "q10": [server.gene_associations("PCSK9", mask="synonymous", maf="<0.1%", limit=1)],
        # Multi-entity questions. The first ten are all single-lookup by
        # construction, which is why they could not see that batch, replication-
        # screen and contrast patterns were missing entirely.
        "q11": [server.catalog("biobanks", trait="T2Diab")],
        # No max_p override: a screen that needs one to show its candidates is
        # the defect, not the call site.
        "q12": [server.phenotype_associations(
            "LDLC", genes="PCSK9,ACAN,TTN", detailed=True, limit=10)],
        "q13": [server.gene_phenotype_detail("PCSK9,LDLR,APOB,ANGPTL3,ABCG5", "LDLC")],
        "q14": [server.top_associations(ancestry="AFR", absent_in="EUR", limit=1)],
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
