"""Run the best optimiser the chains wrote over every beatable Packomania entry, batch by batch.

The chains tuned an optimiser on a 24-case benchmark. This sweep feeds it all the dominated entries
(about 2,000), 60 at a time, on the rented machine, checks every packing it prints in exact
arithmetic, and keeps the best verified packing per entry in runs/packing_best/. No model calls.

    MOONSHOTS_REMOTE=... python3 -m packing.sweep --program runs/packing_best_program.py
"""
import argparse
import json
import os
import subprocess

from common import sandbox
from . import table, verify

ROOT = os.path.dirname(table.ROOT)
BEST = os.path.join(ROOT, "packing_best")
HARVEST = os.path.join(ROOT, "packing_harvest")
REMOTE_DATA = os.environ.get("MOONSHOTS_PACKING_DATA", "/data/packomania")


def all_cases(kinds):
    out = []
    for kind in kinds:
        listed, authors, dom = table.radii(kind), table.authors(kind), table.dominated(kind)
        for n in sorted(dom):
            m, _ = dom[n]
            out.append({"kind": kind, "n": n, "listed_radius": listed[n], "listed_author": authors.get(n, ""),
                        "donor_n": m, "donor_radius": listed[m], "donor_file": f"{REMOTE_DATA}/{kind}_coords/{kind}{m}.txt"})
    return out


def current_best(kind, n):
    for d in (BEST, HARVEST):
        p = os.path.join(d, f"{kind}{n}.json")
        if os.path.exists(p):
            return float(verify.exact(json.load(open(p))["radius"]))
    return 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--program", required=True)
    ap.add_argument("--kinds", nargs="*", default=["cci", "scu", "ssp", "csq"])
    ap.add_argument("--batch", type=int, default=60)
    ap.add_argument("--max-n", type=int, default=3000)
    ap.add_argument("--minutes", type=int, default=30)
    a = ap.parse_args()
    os.makedirs(BEST, exist_ok=True)
    src = open(a.program).read()
    cases = [c for c in all_cases(a.kinds) if c["n"] <= a.max_n]
    done_path = os.path.join(ROOT, "packing_sweep_done.json")
    done = set(json.load(open(done_path))) if os.path.exists(done_path) else set()
    todo = [c for c in cases if f"{c['kind']}{c['n']}" not in done]
    print(f"{len(cases)} cases, {len(todo)} still to sweep", flush=True)
    host, key = os.environ["MOONSHOTS_REMOTE"], os.environ["MOONSHOTS_REMOTE_KEY"]
    for b in range(0, len(todo), a.batch):
        batch = todo[b:b + a.batch]
        path = f"/data/sweep_{b}.json"
        subprocess.run(["ssh", "-i", key, host, f"cat > {path}"], input=json.dumps(batch), text=True, check=True)
        prog = src.replace("/data/packing_cases.json", path).replace("NP=20;", f"NP={min(60, len(batch))};") \
                  .replace("NP=20 ", f"NP={min(60, len(batch))} ").replace("\nNP=20", f"\nNP={min(60, len(batch))}")
        answer, out, err = sandbox.run_full(prog, a.minutes * 60)
        gained = 0
        for p in (answer or {}).get("packings", []) if isinstance(answer, dict) else []:
            try:
                kind, n = str(p["kind"]), int(p["n"])
                if len(p["points"]) != n:
                    continue
                ok, why = verify.verify(kind, p["points"], p["radius"])
            except (KeyError, TypeError, ValueError, ZeroDivisionError):
                continue
            r = float(verify.exact(p["radius"]))
            if ok and r > current_best(kind, n) and r > float(table.radii(kind)[n]):
                listed = table.radii(kind)[n]
                p["listed_radius"], p["gain"] = listed, r / float(listed) - 1
                json.dump(p, open(os.path.join(BEST, f"{kind}{n}.json"), "w"))
                gained += 1
        done.update(f"{c['kind']}{c['n']}" for c in batch)
        json.dump(sorted(done), open(done_path, "w"))
        print(f"batch {b // a.batch + 1}: {len(batch)} cases ({batch[0]['kind']}{batch[0]['n']}..{batch[-1]['kind']}{batch[-1]['n']}), "
              f"{gained} improved beyond the shortcut | {err.strip().splitlines()[-1][:120] if err.strip() else ''}", flush=True)


if __name__ == "__main__":
    main()
