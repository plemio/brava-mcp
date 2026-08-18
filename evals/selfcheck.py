"""Walk each eval question through the tools and check the gold answer comes out.

This is not the full evaluation: it fixes the tool path instead of letting a
model choose it, so it cannot detect a model getting lost. What it does prove is
that every question IS answerable through the exposed surface, and at what call
cost. A server that answers correctly in 25 calls is a bad server, so the call
count is the number to watch, not just the hit rate.

The LLM-in-the-loop runner (thread_ephemeral against the local engine) is the
missing half; see the skill's references/evaluation.md.

    uv run python evals/selfcheck.py
"""

import asyncio
import json
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
        "q09": [server.gene_variants("PCSK9", "LDLC", limit=1)],
        "q10": [server.gene_associations("PCSK9", mask="synonymous", maf="<0.1%", limit=1)],
    }


def found(gold: str, blob: str) -> bool:
    """q05/q07 answer with counts the tools report rather than spell out."""
    if gold in blob:
        return True
    return all(part.strip() in blob for part in gold.split("/"))


async def main() -> int:
    client.reset_counters()
    plans = await paths()
    rows, failures = [], 0

    for q in QUESTIONS:
        before = client.outbound_requests()
        calls = plans[q["id"]]
        blob = "".join(await asyncio.gather(*calls))
        ok = found(q["answer"], blob)
        failures += not ok
        rows.append(
            (q["id"], "PASS" if ok else "FAIL", len(calls),
             client.outbound_requests() - before, len(blob), q["answer"])
        )

    print(f"{'id':5} {'':4} {'calls':>5} {'http':>5} {'chars':>7}  answer")
    for r in rows:
        print(f"{r[0]:5} {r[1]:4} {r[2]:>5} {r[3]:>5} {r[4]:>7}  {r[5]}")
    print()
    print(f"accuracy {len(rows) - failures}/{len(rows)}"
          f" | median calls {sorted(r[2] for r in rows)[len(rows) // 2]}"
          f" | total http {sum(r[3] for r in rows)}"
          f" | total chars {sum(r[4] for r in rows)}")
    await client.close()
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(asyncio.run(main()))
