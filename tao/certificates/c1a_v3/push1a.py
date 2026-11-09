"""C1d upper bound: signed step function v on n equal cells of [-1/4,1/4]; R(v) = 2n max_k |c_k| / S^2, c = v*v, S = sum v.
L-BFGS on a soft maximum of |r_k| (r = 2n c / S^2) with rising sharpness p, re-softened in rounds; FFT for c and gradient.
    python push1d.py in.json out.json [--split K] [--iters N] [--plist ...] [--rounds R]"""
import argparse, json, time
import numpy as np
from scipy.optimize import minimize
from scipy.signal import fftconvolve

def ratio(v):
    c = fftconvolve(v, v); return 2 * len(v) * np.abs(c).max() / v.sum() ** 2

def objective(v, p):
    n = len(v); S = v.sum()
    r = 2 * n * fftconvolve(v, v) / S ** 2
    M = np.abs(r).max()
    ep = np.exp(p * (r - M)); em = np.exp(p * (-r - M)); Z = ep.sum() + em.sum()
    sm = M + np.log(Z) / p
    g = (ep - em) / Z                                   # d sm / d r_k
    corr = fftconvolve(g, v[::-1])[n - 1:2 * n - 1]     # sum_k g_k v_{k-i}
    grad = 2 * n * 2 * corr / S ** 2 - 2 * np.dot(g, r) / S
    return sm, grad

def save(v, path):
    vi = np.rint(np.maximum(v, 0) * (2.0 ** 56 / v.max())).astype(np.int64)
    json.dump({"problem": "C_1a upper bound witness: non-negative step function, n equal cells on [-1/4,1/4], cell values v_int",
               "n": str(len(vi)), "v_int": [int(x) for x in vi]}, open(path, "w"))
    return ratio(vi.astype(float))

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("inp"); ap.add_argument("out")
    ap.add_argument("--split", type=int, default=1); ap.add_argument("--iters", type=int, default=20000)
    ap.add_argument("--plist", default="1e4,1e5,1e6,1e7"); ap.add_argument("--rounds", type=int, default=3)
    a = ap.parse_args()
    d = json.load(open(a.inp)); v = np.array([float(x) for x in d["v_int"]]); v = np.repeat(v, a.split)
    v /= np.abs(v).max()
    best_v, best = v.copy(), ratio(v); print(f"n={len(v)} start={best:.13f}", flush=True); t0 = time.time()
    for rnd in range(a.rounds):
        for p in [float(x) for x in a.plist.split(",")]:
            res = minimize(objective, best_v, args=(p,), jac=True, method="L-BFGS-B", bounds=[(0, None)] * len(best_v),
                           options={"maxiter": a.iters, "maxcor": 50, "ftol": 0, "gtol": 0})
            v = res.x / np.abs(res.x).max(); s = ratio(v)
            si = float("nan")
            if s < best:
                best_v, best = v.copy(), s; si = save(best_v, a.out)
            print(f"r{rnd} p={p:.0e} it={res.nit} ratio={s:.13f} best={best:.13f} int={si:.13f} t={time.time()-t0:.0f}s", flush=True)

if __name__ == "__main__":
    main()
