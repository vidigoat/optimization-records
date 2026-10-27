"""LLL generator v2: clustered dual measure (bump centres + touch points), exact high-precision evaluation lattice."""
import numpy as np, mpmath as mp
from flint import fmpz_mat, fmpz_poly

def clusters(t, nu, g):
    """Aggregate dual measure into atoms: contiguous runs of support -> (t_centre (mass-weighted), mass)."""
    sup = np.flatnonzero(nu > 0)
    groups = np.split(sup, np.flatnonzero(np.diff(sup) > 1) + 1)
    out = []
    for gr in groups:
        m = nu[gr].sum(); tc = (nu[gr] @ t[gr]) / m
        out.append((tc, m, t[gr[-1]] - t[gr[0]]))
    return np.array(out)   # cols: t, mass, width

def lattice_polys(ts, wts, d, lam=1.0, bits=200, nkeep=None, profile=None):
    """Short integer vectors for  lam^2 |c|^2 + sum wts |S(w)|^2 * exp(-2 d profile)."""
    mp.mp.prec = bits + 64
    S = mp.mpf(2) ** bits
    rows = [[0] * (d + 1) for _ in range(d + 1)]
    cols = []
    for j, (tt, ww) in enumerate(zip(ts, wts)):
        z = mp.expjpi(2 * mp.mpf(tt)); w = z * (1 - z)
        sc = mp.sqrt(mp.mpf(ww)) * (mp.exp(-d * mp.mpf(profile[j])) if profile is not None else 1)
        p = sc; re = []; im = []
        for k in range(d + 1):
            re.append(int(mp.nint(p.real * S))); im.append(int(mp.nint(p.imag * S)))
            p = p * w
        cols.append(re); cols.append(im)
    lamI = int(mp.nint(mp.mpf(lam) * S))
    M = []
    for k in range(d + 1):
        M.append([lamI if i == k else 0 for i in range(d + 1)] + [c[k] for c in cols])
    L = fmpz_mat(M).lll()
    out = []
    nk = nkeep or (d + 1)
    for r in range(nk):
        v = [int(L[r, k]) for k in range(d + 1)]
        q = [x // lamI if x % lamI == 0 else None for x in v]
        if None in q: continue
        if not any(q[1:]): continue
        out.append(q)       # ascending coefficients
    return out, L

def factors(coeffs_asc):
    P = fmpz_poly(coeffs_asc)
    res = []
    for f, e in P.factor()[1]:
        if f.degree() >= 1:
            c = [int(x) for x in f.coeffs()][::-1]
            if c[0] < 0: c = [-x for x in c]
            res.append(c)
    return res
