"""Export a C86 search state to a certificate in R. Pinho's format: exact rational edges and node densities
(piecewise linear, zero outside), auxiliary polynomials A (constant first) with rational weights b, bound = T/m.

Knots are the exact rationals of the binary floats used in the search; triangle weights are rounded to multiples
of 1e-15 and auxiliary weights rounded down to multiples of 1e-14. Only the frozen certificate is a claim; it must
pass the exact checkers.
    python export86.py state.pkl out.json
"""
import sys, json, pickle
from fractions import Fraction as F
import numpy as np


def export(st, out):
    tris, w = st['tris'], st['w']
    keep = w > 1e-13
    T, W = tris[keep], w[keep]
    knots = sorted({F(float(v)) for v in T.ravel()})
    pos = {k: i for i, k in enumerate(knots)}
    dens = [F(0)] * len(knots)
    for (l, c, r), wt in zip(T, W):
        l, c, r = F(float(l)), F(float(c)), F(float(r))
        wq = F(round(float(wt) * 10**15), 10**15)
        if wq == 0: continue
        peak = 2 * wq / (r - l)                               # unit-mass triangle scaled by wq
        il, ic, ir = pos[l], pos[c], pos[r]
        for i in range(il + 1, ir):
            k = knots[i]
            dens[i] += peak * ((k - l) / (c - l) if k <= c else (r - k) / (r - c))
    # zero-density endpoints are kept so the density is continuous
    m = sum((knots[i + 1] - knots[i]) * (dens[i] + dens[i + 1]) / 2 for i in range(len(knots) - 1))
    Tm = sum((knots[i + 1] - knots[i]) * (dens[i] * (2 * knots[i] + knots[i + 1]) + dens[i + 1] * (knots[i] + 2 * knots[i + 1])) / 6
             for i in range(len(knots) - 1))
    A, B = [], []
    merged = {}
    for Q, bj in zip(st['polys'], st['b']):
        Q = [int(t) for t in Q]
        if Q[-1] < 0: Q = [-t for t in Q]
        merged[tuple(Q)] = merged.get(tuple(Q), 0.0) + max(float(bj), 0.0)   # identical columns: weights add
    for Q, bj in merged.items():
        bq = F(int(bj * 10**14), 10**14)
        if bq > 0: A.append(list(Q)); B.append(bq)
    degsum = sum(b * (len(a) - 1) for a, b in zip(A, B))
    assert degsum < m, 'degree condition'
    bound = Tm / m
    cert = {'edges': [str(k) for k in knots], 'node_densities': [str(d) for d in dens], 'A': A, 'b': [str(b) for b in B],
            'bound': str(bound), 'bound_float': float(bound)}
    json.dump(cert, open(out, 'w'))
    print(f'{len(knots)} knots, {len(A)} auxiliaries, mass {float(m):.15f}, sum b deg / m {float(degsum / m):.6f}, bound {float(bound):.12f}')
    return cert


if __name__ == '__main__':
    st = pickle.load(open(sys.argv[1], 'rb'))
    export(st, sys.argv[2])
