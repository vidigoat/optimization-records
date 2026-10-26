"""Discretised Talebizadeh Sardari–Orloski primal LP for C86 (Schur–Siegel–Smyth trace constant), as in
Koshta's c86-sss-trace-certificate, with any list of auxiliary polynomials.

Variables: densities d_i >= 0 on cells [e_i, e_{i+1}] inside [a, b], weights b_j >= 0.
Minimise T = sum d_i (e_{i+1}^2 - e_i^2)/2 subject to
  m = sum d_i (e_{i+1} - e_i) = 1,
  f(x_s) = U_mu(x_s) - sum_j b_j log|Q_j(x_s)| >= margin at sample points x_s in [a, b],
  sum_j b_j deg Q_j <= 1 - slack,
  int log|Q_j| dmu >= 0 for every j.
Writes a certificate in the [Pin26]/[Kos2026] JSON format (rationals), to be checked by checkers/c86_trace.py.

    python lp.py polys.txt out.json [--h 200] [--a 11/200] [--b 519/100] [--sub 4] [--margin 2e-7]
"""
import argparse
import ast
import json
from fractions import Fraction

import numpy as np
from scipy.optimize import linprog


def G(u):
    u = np.asarray(u, dtype=float)
    out = np.zeros_like(u)
    nz = u != 0
    out[nz] = u[nz] * np.log(np.abs(u[nz])) - u[nz]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("polys")
    ap.add_argument("out")
    ap.add_argument("--h", type=int, default=200)          # cells of width 1/h
    ap.add_argument("--a", default="11/200")
    ap.add_argument("--b", default="519/100")
    ap.add_argument("--sub", type=int, default=4)          # sample points per cell
    ap.add_argument("--margin", type=float, default=2e-7)
    ap.add_argument("--slack", type=float, default=1e-9)
    ap.add_argument("--maxdeg", type=int, default=99)
    ap.add_argument("--e4", type=float, default=1e-10)
    ap.add_argument("--prune", type=int, default=0)
    a = ap.parse_args()
    A, B = Fraction(a.a), Fraction(a.b)
    edges = []
    k = 0
    while True:
        e = A + Fraction(k, a.h)
        if e > B:
            break
        edges.append(e)
        k += 1
    if edges[-1] != B:
        edges.append(B)
    E = np.array([float(e) for e in edges])
    N = len(E) - 1
    polys = []
    for line in open(a.polys):
        line = line.strip()
        if line and not line.startswith("#"):
            c = ast.literal_eval(line)
            if len(c) - 1 <= a.maxdeg:
                polys.append(c)
    P = len(polys)
    print(f"cells={N} polys={P}", flush=True)
    roots = [np.roots(c[::-1]).real for c in polys]
    degs = np.array([len(c) - 1 for c in polys])
    # sample points
    hcell = 1.0 / a.h
    near = hcell * np.array([1e-6, 1e-5, 1e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.1])
    pts = [np.linspace(E[i], E[i + 1], a.sub, endpoint=False) for i in range(N)] + [[E[-1]]]
    pts += [E[:, None] + near[None, :], E[:, None] - near[None, :]]
    xs = np.unique(np.concatenate([np.ravel(p) for p in pts]))
    xs = xs[(xs >= E[0]) & (xs <= E[-1])]
    ILQ = np.zeros((P, N))
    for j, rs in enumerate(roots):
        for r in rs:
            ILQ[j] += G(E[1:] - r) - G(E[:-1] - r)
    widths = E[1:] - E[:-1]
    moments = (E[1:] ** 2 - E[:-1] ** 2) / 2
    c = np.concatenate([moments, np.zeros(P)])

    def rows(x):
        Um = G(x[:, None] - E[None, :-1]) - G(x[:, None] - E[None, 1:])
        lq = np.zeros((len(x), P))
        for j, cf in enumerate(polys):
            lq[:, j] = np.log(np.maximum(np.abs(np.polyval(cf[::-1], x)), 1e-300))
        return np.hstack([-Um, lq])

    if a.prune:
        # coarse solve (2 samples per cell) to rank polynomials by reduced cost; keep the best a.prune
        xc = np.concatenate([np.linspace(E[i], E[i + 1], 2, endpoint=False) for i in range(N)] + [[E[-1]]])
        A_ub = np.vstack([rows(xc), np.hstack([np.zeros(N), degs])[None, :], np.hstack([-ILQ, np.zeros((P, P))])])
        b_ub = np.concatenate([np.full(len(xc), -a.margin), [1 - a.slack], np.full(P, -a.e4)])
        A_eq = np.hstack([widths, np.zeros(P)])[None, :]
        rc = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=[1.0], bounds=[(0, None)] * (N + P), method="highs")
        red = rc.lower.marginals[N:]
        keep = np.argsort(red)[:a.prune]
        keep = np.union1d(keep, np.where(rc.x[N:] > 1e-12)[0])
        print(f"coarse LP {rc.fun:.10f}; keeping {len(keep)} of {P} polynomials", flush=True)
        polys = [polys[j] for j in keep]
        roots = [roots[j] for j in keep]
        degs = degs[keep]
        ILQ = ILQ[keep]
        P = len(polys)
        c = np.concatenate([moments, np.zeros(P)])
    fine = np.concatenate([np.linspace(E[i], E[i + 1], 200, endpoint=False) for i in range(N)] + [[E[-1]]])
    fine = np.unique(np.concatenate([fine, xs]))
    Fine = None
    for rnd in range(30):
        A_ub = np.vstack([rows(xs),
                          np.hstack([np.zeros(N), degs])[None, :],
                          np.hstack([-ILQ, np.zeros((P, P))])])
        b_ub = np.concatenate([np.full(len(xs), -a.margin), [1 - a.slack], np.full(P, -a.e4)])
        A_eq = np.hstack([widths, np.zeros(P)])[None, :]
        res = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=[1.0], bounds=[(0, None)] * (N + P), method="highs")
        # f on the fine grid, in chunks
        fv = np.concatenate([-(rows(fine[k:k + 20000]) @ res.x) for k in range(0, len(fine), 20000)])
        bad = fine[fv < a.margin / 2]
        print(f"round {rnd}: LP {res.fun:.10f}  samples {len(xs)}  fine-grid min f {fv.min():.3e}  violations {len(bad)}", flush=True)
        if len(bad) == 0:
            break
        xs = np.unique(np.concatenate([xs, bad]))
    d, bw = res.x[:N], res.x[N:]
    print(f"LP bound T/m = {res.fun:.10f}; polys used: {(bw > 1e-12).sum()}", flush=True)
    # rational rounding: densities to 1e-15 grid, weights rounded down to 1e-15
    dq = [Fraction(int(round(x * 1e15)), 10 ** 15) if x > 1e-14 else Fraction(0) for x in d]
    bq = [Fraction(int(np.floor(x * 1e15)), 10 ** 15) if x > 1e-14 else Fraction(0) for x in bw]
    used = [j for j in range(P) if bq[j] > 0]
    m = sum(dq[i] * (edges[i + 1] - edges[i]) for i in range(N))
    T = sum(dq[i] * (edges[i + 1] ** 2 - edges[i] ** 2) / 2 for i in range(N))
    cert = {"edges": [str(e) for e in edges], "densities": [str(x) for x in dq],
            "A": [polys[j] for j in used], "b": [str(bq[j]) for j in used], "bound": str(T / m),
            "bound_float": float(T / m)}
    json.dump(cert, open(a.out, "w"))
    print(f"wrote {a.out}: bound {float(T / m):.12f} with {len(used)} polynomials", flush=True)
    for j in used:
        if len(polys[j]) - 1 >= 7:
            print(f"  deg {len(polys[j]) - 1}: {polys[j]}  b={float(bq[j]):.3e}", flush=True)


if __name__ == "__main__":
    main()
