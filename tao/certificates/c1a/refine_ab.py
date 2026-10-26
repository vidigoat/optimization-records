"""Refine step functions for C1a (Sidon autocorrelation) and C1b (Erdős minimum overlap), Tao's table, upper bounds.

1a: heights a >= 0, n cells on [-1/4, 1/4];  score = 2n max_k (a*a)_k / (sum a)^2        (minimise)
1b: heights 0 <= h <= 1, n cells on [0, 2], sum h = n/2;  score = 2 max_s c_s / n,
    c = h * reverse(1 - h)  (full convolution, length 2n-1)                                  (minimise)
L-BFGS-B on a soft maximum with p increasing; FFT for the correlations and gradients. For 1b the mass is held
by a quadratic penalty and restored exactly at the end.

    python refine_ab.py 1a|1b in.json split out.json [--plist ...] [--iters N] [--check]
"""
import argparse
import json
from fractions import Fraction

import numpy as np
from scipy.optimize import minimize
from scipy.signal import fftconvolve


def load(path):
    vals = json.load(open(path))["values"]
    return np.array([float(Fraction(v)) for v in vals])


def score(mode, x):
    n = len(x)
    if mode == "1a":
        return 2 * n * fftconvolve(x, x).max() / x.sum() ** 2
    c = fftconvolve(x, (1 - x)[::-1])
    return 2 * c.max() / n


def objective(x, mode, p, lam):
    n = len(x)
    if mode == "1a":
        S = x.sum()
        c = fftconvolve(x, x)
        r = 2 * n * c / S ** 2
        M = r.max()
        e = np.exp(p * (r - M))
        Z = e.sum()
        sm = M + np.log(Z) / p
        w = e / Z
        # d r_m / d x_i = 2n * 2 x_{m-i} / S^2 - 2 r_m / S
        g = 2 * n * 2 * fftconvolve(w, x[::-1], mode="valid") / S ** 2 - 2 * np.dot(w, r) / S
        return sm, g
    g1 = 1 - x
    c = fftconvolve(x, g1[::-1])                 # c_s, s = 0..2n-2 ; c_s = sum_i x_i g1_{i+n-1-s}
    r = 2 * c / n
    M = r.max()
    e = np.exp(p * (r - M))
    Z = e.sum()
    sm = M + np.log(Z) / p
    w = e / Z * 2 / n                            # d sm / d c_s
    # c_s = sum_i x_i (1 - x_{i+n-1-s});  d/dx_k: w_s (1 - x_{k+n-1-s}) - w_s x_{k-(n-1-s)}
    t1 = fftconvolve(w[::-1], g1, mode="full")   # sum_s w_s g1_{k+n-1-s}
    t1 = t1[n - 1:2 * n - 1] if False else _corr_g(w, g1, n)
    t2 = _corr_h(w, x, n)
    grad = t1 - t2
    pen = lam * (x.sum() - n / 2) ** 2
    return sm + pen, grad + 2 * lam * (x.sum() - n / 2)


def _corr_g(w, g1, n):
    # T1_k = sum_s w_s g1_{k + n - 1 - s}; let u = n-1-s (u from -(n-1) .. n-1), w'_u = w_{n-1-u}
    # T1_k = sum_u w'_u g1_{k+u}
    full = fftconvolve(g1, w, mode="full")       # (g1*w)_t = sum_s g1_{t-s} w_s ; set t = k + n - 1
    return full[n - 1:2 * n - 1]


def _corr_h(w, x, n):
    # T2_k = sum_s w_s x_{k - (n-1-s)} = sum_s w_s x_{k - n + 1 + s}
    full = fftconvolve(x, w[::-1], mode="full")  # sum_s x_{t - (2n-2-s')}... index via valid range
    return full[n - 1:2 * n - 1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode")
    ap.add_argument("inp")
    ap.add_argument("split", type=int)
    ap.add_argument("out")
    ap.add_argument("--plist", default="1e5,3e5,1e6,3e6,1e7,3e7")
    ap.add_argument("--iters", type=int, default=3000)
    ap.add_argument("--lam", type=float, default=1e3)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    x = load(a.inp)
    print(f"input n={len(x)} score={score(a.mode, x):.12f}", flush=True)
    x = np.repeat(x, a.split)
    if a.mode == "1a":
        x = x / x.sum()
    print(f"split n={len(x)} score={score(a.mode, x):.12f}", flush=True)
    if a.check:
        rng = np.random.default_rng(0)
        y = rng.random(50) * 0.8 + 0.1
        for i in [0, 13, 49]:
            e = np.zeros(50); e[i] = 1e-6
            f0, g = objective(y, a.mode, 30.0, 1.0)
            fd = (objective(y + e, a.mode, 30.0, 1.0)[0] - objective(y - e, a.mode, 30.0, 1.0)[0]) / 2e-6
            print(f"grad check {i}: {g[i]:.8f} vs fd {fd:.8f}")
        return
    best, bx = score(a.mode, x), x.copy()
    bounds = [(0, None)] * len(x) if a.mode == "1a" else [(0, 1)] * len(x)
    for p in [float(v) for v in a.plist.split(",")]:
        res = minimize(objective, x, args=(a.mode, p, a.lam), jac=True, method="L-BFGS-B", bounds=bounds,
                       options={"maxiter": a.iters, "maxcor": 30, "ftol": 0, "gtol": 0})
        x = res.x
        if a.mode == "1a":
            x = x / x.sum()
        s = score(a.mode, x)
        if a.mode == "1b":
            s_mass = x.sum() - len(x) / 2
        else:
            s_mass = 0.0
        if s < best:
            best, bx = s, x.copy()
        print(f"  p={p:.0e} score={s:.12f} best={best:.12f} mass_err={s_mass:.2e}", flush=True)
    json.dump({"values": [repr(float(v)) for v in bx]}, open(a.out, "w"))
    print(f"saved {a.out}", flush=True)


if __name__ == "__main__":
    main()
