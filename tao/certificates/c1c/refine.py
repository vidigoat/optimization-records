"""Refine a step function for C1c (Barnard–Steinerberger autocorrelation constant, Tao's table, lower bound).

Score of non-negative weights w on cells of width 1/K:  K * min_{0<=k<=K} A_k / (sum w)^2,  A_k = sum_i w_i w_{i+k}.
Maximised with L-BFGS-B on a soft minimum, p increasing; FFT for A and for the gradient.

    python refine.py in.json split out.json [--iters N]
"""
import argparse
import json

import numpy as np
from scipy.optimize import minimize
from scipy.signal import fftconvolve


def autocorr(w):
    return fftconvolve(w, w[::-1])[len(w) - 1:]          # A_k for k = 0..n-1


def score(w, K):
    A = autocorr(w)[:K + 1]
    return K * A.min() / w.sum() ** 2


def objective(w, K, p):
    n = len(w)
    S = w.sum()
    A = autocorr(w)[:K + 1]
    r = K * A / S ** 2
    rmin = r.min()
    e = np.exp(-p * (r - rmin))
    Z = e.sum()
    sm = rmin - np.log(Z / len(r)) / p                    # soft minimum, <= true min + log(len)/p
    g = e / Z                                              # d sm / d r_k
    # d r_k / d w_i = K (w_{i+k} + w_{i-k}) / S^2 - 2 r_k / S
    gk = np.zeros(n)
    gk[:K + 1] = g
    t1 = fftconvolve(w, gk[::-1])[n - 1:2 * n - 1]         # sum_k g_k w_{i+k}
    t2 = fftconvolve(w, gk)[:n]                            # sum_k g_k w_{i-k}
    grad = K * (t1 + t2) / S ** 2 - 2 * np.dot(g, r) / S
    return -sm, -grad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inp")
    ap.add_argument("split", type=int)
    ap.add_argument("out")
    ap.add_argument("--iters", type=int, default=3000)
    ap.add_argument("--plist", default="2e5,5e5,1e6,2e6,5e6,1e7")
    a = ap.parse_args()
    d = json.load(open(a.inp))
    K0 = d["K"]
    w = np.array(d["weights_int"], dtype=float)
    w /= w.sum()
    print(f"input K={K0} n={len(w)} score={score(w, K0):.12f}", flush=True)
    K = K0 * a.split
    w = np.repeat(w, a.split) / a.split
    print(f"split K={K} n={len(w)} score={score(w, K):.12f}", flush=True)
    best_w, best = w.copy(), score(w, K)
    for p in [float(x) for x in a.plist.split(",")]:
        res = minimize(objective, w, args=(K, p), jac=True, method="L-BFGS-B",
                       bounds=[(0, None)] * len(w), options={"maxiter": a.iters, "maxcor": 30, "ftol": 0, "gtol": 0})
        w = res.x / res.x.sum()
        s = score(w, K)
        if s > best:
            best_w, best = w.copy(), s
        print(f"  p={p:.0e} score={s:.12f} best={best:.12f}", flush=True)
    scale = 2.0 ** 40 / best_w.max()
    ints = np.rint(best_w * scale).astype(np.int64)
    json.dump({"K": K, "weights_int": [int(x) for x in ints]}, open(a.out, "w"))
    print(f"saved {a.out}: float score of integer weights {score(ints.astype(float), K):.12f}", flush=True)


if __name__ == "__main__":
    main()
