"""L-BFGS refinement of a symmetric step graphon for the K4 Ramsey multiplicity (Tao's table, C89).

Same parametrisation as JeetuSK0808's refine.py: blocks grouped into classes by (sorted row, diagonal), every
off-diagonal pair keyed by (class i, class j, original value); pairs with one key move together, so the
density and gradient come from one representative block per class. Diagonal entries are fixed.
Search only (float64): the output is rounded to multiples of 2^-DEN_BITS and must be certified by exact_eval.py.

usage: refine_lbfgs.py in.json.gz out.json.gz iterations [--check]
"""
import gzip
import json
import sys
import time

import numpy as np
from scipy.optimize import minimize

d = json.load(gzip.open(sys.argv[1]))
D = d["edge_probability_denominator"]
R0 = np.array(d["red_probability_numerators"], dtype=np.int64)
m = R0.shape[0]
srt = np.sort(R0, axis=1)
_, cls = np.unique(np.hstack([srt, np.diag(R0)[:, None]]), axis=0, return_inverse=True)
cls = cls.ravel()
cid = range(int(cls.max()) + 1)
nc = len(cid)
reps = np.array([np.nonzero(cls == a)[0][0] for a in range(nc)])
size = np.bincount(cls)
a_ = np.minimum(cls[:, None], cls[None, :]); b_ = np.maximum(cls[:, None], cls[None, :])
code = (a_ * nc + b_) * (int(R0.max()) + 1) + R0
code[np.arange(m), np.arange(m)] = -1
uniq, inv = np.unique(code, return_inverse=True)
inv = inv.reshape(m, m)
K = np.where(code >= 0, inv - (1 if uniq[0] == -1 else 0), -1)
keyid = range(int(K.max()) + 1)
nk = len(keyid)
off = K >= 0
Kc = np.where(off, K, 0)
x0 = np.zeros(nk)
x0[K[off]] = (R0 / D)[off]
diagW = np.diag(R0 / D).copy()
print(f"m={m} classes={nc} keys={nk}", flush=True)


def build(x):
    W = x[Kc]
    np.fill_diagonal(W, diagW)
    return W


def density_grad(x):
    W = build(x)
    tot = 0.0
    grow = np.zeros((nc, m))
    for col, sgn in ((W, 1.0), (1.0 - W, -1.0)):
        for a in range(nc):
            i = reps[a]
            A = col * col[i][None, :]
            Q = A @ col
            S = np.sum(Q * A, axis=1)
            tot += size[a] * np.dot(col[i], S)
            grow[a] += sgn * S
    g = np.zeros(nk)
    for a in range(nc):
        i = reps[a]
        k = K[i]
        ok = k >= 0
        np.add.at(g, k[ok], size[a] * grow[a][ok])
    return tot / m ** 4, g / 2.0 * 12.0 / m ** 4


c0, g0 = density_grad(x0)
print(f"start c4={c0:.15f}", flush=True)

if "--check" in sys.argv:          # finite-difference check of the gradient on a few keys
    rng = np.random.default_rng(1)
    for k in rng.choice(nk, 3, replace=False):
        h = 1e-6
        xp = x0.copy(); xp[k] += h
        xm = x0.copy(); xm[k] -= h
        fd = (density_grad(xp)[0] - density_grad(xm)[0]) / (2 * h)
        print(f"key {k}: analytic {g0[k]:.6e}  finite-diff {fd:.6e}", flush=True)
    sys.exit()

SCALE = 1e9
best = [c0, x0.copy()]
t0 = time.time()


def fun(x):
    c, g = density_grad(x)
    if c < best[0]:
        best[0], best[1] = c, x.copy()
        best.append(1)
        if len(best) % 25 == 0:
            save(x)
    print(f"  eval c4={c:.15f} gain={c0 - c:.3e} ({time.time() - t0:.0f}s)", flush=True)
    return (c - c0) * SCALE, g * SCALE


def save(x, bits=30):
    DO = 1 << bits
    W = build(x)
    Rr = np.rint(W * DO).astype(np.int64)
    Rr = np.minimum(Rr, Rr.T)
    out = dict(d)
    out["edge_probability_denominator"] = DO
    out["red_probability_numerators"] = Rr.tolist()
    json.dump(out, gzip.open(sys.argv[2], "wt"))


res = minimize(fun, x0, jac=True, method="L-BFGS-B", bounds=[(0.0, 1.0)] * nk,
               options={"maxiter": int(sys.argv[3]), "maxfun": int(sys.argv[3]) * 2, "ftol": 0, "gtol": 0, "maxcor": 30})
save(best[1])
print(f"done: c4={best[0]:.15f} gain={c0 - best[0]:.3e} ({res.message})", flush=True)
