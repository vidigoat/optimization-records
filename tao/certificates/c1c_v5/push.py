"""Push a C1c step function: optional cell split, then L-BFGS-B on a soft minimum with rising sharpness p.
Saves the best integer certificate after every stage.

    python push.py in.json split out.json [--iters N] [--plist ...]
"""
import argparse, json, time
import numpy as np
from scipy.optimize import minimize
from refine import score, objective

def save(w, K, path):
    ints = np.rint(w * (2.0 ** 56 / w.max())).astype(np.int64)
    json.dump({"K": K, "weights_int": [int(x) for x in ints]}, open(path, "w"))
    return score(ints.astype(float), K)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inp"); ap.add_argument("split", type=int); ap.add_argument("out")
    ap.add_argument("--iters", type=int, default=5000)
    ap.add_argument("--plist", default="3e7,1e8,3e8,1e9")
    ap.add_argument("--rounds", type=int, default=1)
    a = ap.parse_args()
    d = json.load(open(a.inp)); K = d["K"] * a.split
    w = np.repeat(np.array(d["weights_int"], dtype=float), a.split); w /= w.sum()
    best_w, best = w.copy(), score(w, K)
    print(f"K={K} n={len(w)} start={best:.13f}", flush=True); t0 = time.time()
    for rnd in range(a.rounds):
        for p in [float(x) for x in a.plist.split(",")]:
            res = minimize(objective, best_w, args=(K, p), jac=True, method="L-BFGS-B", bounds=[(0, None)] * len(w),
                           options={"maxiter": a.iters, "maxcor": 50, "ftol": 0, "gtol": 0})
            w = res.x / res.x.sum(); s = score(w, K)
            if s > best:
                best_w, best = w.copy(), s
                si = save(best_w, K, a.out)
            else:
                si = float("nan")
            print(f"r{rnd} p={p:.0e} it={res.nit} score={s:.13f} best={best:.13f} int={si:.13f} t={time.time()-t0:.0f}s", flush=True)

if __name__ == "__main__":
    main()
