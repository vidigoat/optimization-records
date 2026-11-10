"""Search for new auxiliary polynomials for the C86 LP by iterated dual-weighted lattice reduction.

A new column b_Q >= 0 lowers the LP value only if its reduced cost is negative:
    rc(Q) = sum_i y_i log|Q(x_i)| + z deg Q  < 0,
with y_i >= 0 the duals of the potential rows at the sample points x_i and z >= 0 the dual of the degree row.
Minimising sum_i W_i |Q(x_i)|^2 over integer coefficient vectors is a lattice problem; with W_i = y_i / |Q_0(x_i)|^2
it is a quadratic model of the log objective around Q_0, so the step is iterated (IRLS), as in the C82 search.
"""
import numpy as np
from flint import fmpz_mat, fmpz_poly


def rc(cs_low, x, y, z):
    v = np.abs(np.polynomial.polynomial.polyval(x, np.asarray(cs_low, float)))
    if np.any(v == 0):
        return np.inf
    return float(y @ np.log(v) + z * (len(cs_low) - 1))


def lattice_vectors(x, W, d, lam, bits=110):
    """LLL-reduce rows  b_k = (lamS e_k) || (round(S sqrt(W_j) x_j^k))_j,  k = 0..d  (as in the C82 search).
    sum_k c_k b_k has squared length S^2 (lam^2 |c|^2 + sum_j W_j Q(x_j)^2), Q = sum_k c_k x^k."""
    S = 2.0 ** bits
    sw = np.sqrt(W)
    lamI = max(1, int(round(lam * S)))
    M = []
    for k in range(d + 1):
        M.append([lamI if i == k else 0 for i in range(d + 1)] + [int(round(v)) for v in S * sw * x ** k])
    try:
        R = fmpz_mat(M).lll()
    except Exception:
        return []
    out = []
    for r in range(d + 1):
        v = [int(R[r, k]) for k in range(d + 1)]
        if any(t % lamI for t in v) or not any(v): continue
        out.append([t // lamI for t in v])
    return out


def irr_factors(cs_low):
    out = []
    P = fmpz_poly(cs_low)
    if P.degree() < 1:
        return out
    for f, _ in P.factor()[1]:
        c = [int(t) for t in f.coeffs()]
        if c[-1] < 0: c = [-t for t in c]
        if len(c) >= 2: out.append(c)
    return out


def search(x, y, z, starts, rng, lib, libw, dmin=8, dmax=40, iters=15, nk=8):
    """Iterated dual-weighted LLL from products of library polynomials (constant-first lists, weights libw).
    Returns candidates (rc, coeffs) with negative reduced cost on the given duals."""
    found = {}
    p = np.asarray(libw, float) + 1e-12; p = 0.6 * p / p.sum() + 0.4 / len(lib)
    for s in range(starts):
        D0 = int(rng.integers(dmin, dmax + 1))
        Q0 = np.array([1.0]); cs0 = [1]
        for _ in range(30):
            f = lib[int(rng.choice(len(lib), p=p))]
            if len(cs0) - 1 + len(f) - 1 <= D0:
                cs0 = [int(t) for t in (fmpz_poly(cs0) * fmpz_poly(f)).coeffs()]
        if len(cs0) < 4: continue
        Q0 = np.array(cs0, float)
        pert = np.exp(rng.uniform(0, 0.6) * rng.standard_normal(len(x))) if rng.random() < 0.6 else np.ones(len(x))
        extra = int(rng.integers(0, 6)); lam_rel = 10 ** rng.uniform(-6, -1)
        best = rc(cs0, x, y, z)
        for it in range(iters):
            d = len(Q0) - 1 + extra
            v = np.abs(np.polynomial.polynomial.polyval(x, Q0))
            W = y * pert / np.maximum(v, 1e-9 * np.median(v)) ** 2
            lam = lam_rel * np.sqrt(y.sum()) / np.linalg.norm(Q0)
            B = lattice_vectors(x, W, d, lam)
            cands = []
            for i, b in enumerate(B[:nk]):
                cands.append(b)
                for b2 in B[:i]:
                    cands.append([u + w for u, w in zip(b, b2)]); cands.append([u - w for u, w in zip(b, b2)])
            step = None
            for q in cands:
                for f in irr_factors(q):
                    r = rc(f, x, y, z)
                    if not np.isfinite(r): continue
                    if r < -1e-12: found[tuple(f)] = r
                c2 = [int(t) for t in q]
                while c2 and c2[-1] == 0: c2.pop()
                if len(c2) < 2: continue
                if c2[-1] < 0: c2 = [-t for t in c2]
                r2 = rc(c2, x, y, z)
                if np.isfinite(r2) and (step is None or r2 < step[0]) and c2 != cs0: step = (r2, c2)
            if step is None or step[0] >= best - 1e-13: break
            best, cs0 = step; Q0 = np.array(cs0, float)
    return sorted(((r, list(f)) for f, r in found.items()), key=lambda t: t[0])
