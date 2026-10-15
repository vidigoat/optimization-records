"""Turn every record that a larger packing beats into a verified packing of our own.

For each dominated entry n with donor m: take the donor's published centres, delete m - n of them
(the ones whose removal frees the most room is a refinement for the optimiser; here we delete those
nearest the most crowded contacts), set the radius just below what the remaining centres allow,
and check the result in exact arithmetic. Only gains above MIN_GAIN (relative) are kept, so that
nothing rests on rounding. The optimiser chains then improve these further.

    python3 -m packing.harvest --kinds csq cci scu ssp
"""
import argparse
import json
import math
import os
import sys
import time

from . import table, verify

OUT = os.path.join(os.path.dirname(table.ROOT), "packing_harvest")
MIN_GAIN = 1e-9


def tightest(kind, pts):
    """Index of the object in the tightest spot (smallest gap to a neighbour or the wall)."""
    fl = [tuple(float(c) for c in p) for p in pts]
    best, worst = None, math.inf
    r = verify.largest_radius(kind, pts)
    for i, p in enumerate(fl):
        if kind in ("csq", "scu"):
            g = min(0.5 - abs(c) for c in p) - r
        else:
            g = 1 - math.sqrt(sum(c * c for c in p)) - r
        if g < worst:
            worst, best = g, i
    return best


def improve(kind, n, m):
    pts = table.coords(kind, m)
    if len(pts) != m:
        return None, f"donor file has {len(pts)} points, expected {m}"
    while len(pts) > n:                       # remove the object pinned hardest against the wall
        del pts[tightest(kind, pts)]
    r = verify.largest_radius(kind, pts) * (1 - 1e-12)
    rs = repr(r)
    ok, why = verify.verify(kind, pts, rs)
    if not ok:
        return None, why
    return {"kind": kind, "n": n, "radius": rs, "points": [list(p) for p in pts], "donor": m}, why


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kinds", nargs="*", default=["csq", "cci", "scu", "ssp"])
    ap.add_argument("--max-n", type=int, default=100000)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    for kind in a.kinds:
        listed = table.radii(kind)
        dom = table.dominated(kind)
        done = kept = 0
        t0 = time.time()
        for n in sorted(dom):
            if n > a.max_n:
                break
            path = os.path.join(OUT, f"{kind}{n}.json")
            if os.path.exists(path):
                done += 1; kept += 1
                continue
            m, _ = dom[n]
            res, why = improve(kind, n, m)
            done += 1
            if res is None:
                print(f"{kind}{n}: failed ({why})", flush=True)
                continue
            gain = float(verify.exact(res["radius"]) / verify.exact(listed[n])) - 1
            if gain <= MIN_GAIN:
                continue
            res["listed_radius"], res["gain"] = listed[n], gain
            with open(path, "w") as f:
                json.dump(res, f)
            kept += 1
            if kept % 50 == 0:
                print(f"{kind}: {kept} verified improvements so far ({done} tried, {time.time() - t0:.0f}s)", flush=True)
        print(f"{kind}: done, {kept} verified improvements out of {len(dom)} dominated entries", flush=True)


if __name__ == "__main__":
    main()
