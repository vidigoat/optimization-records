"""IRLS lattice local search for integer S minimising rc(S) = sum nu log|S(w)| over the dual support."""
import numpy as np, mpmath as mp
from flint import fmpz_mat, fmpz_poly

class Support:
    def __init__(self, t, nu, sub=3):
        sup = np.flatnonzero(nu > 0)
        groups = np.split(sup, np.flatnonzero(np.diff(sup) > 1) + 1)
        pts, wts = [], []
        for gr in groups:
            if len(gr) <= sub:
                pts += list(t[gr]); wts += list(nu[gr])
            else:
                # sub points evenly inside the run, each carrying 1/sub of the mass
                idx = gr[np.linspace(0, len(gr)-1, sub).round().astype(int)]
                m = nu[gr].sum()
                pts += list(t[idx]); wts += [m/sub]*sub
        self.t = np.array(pts); self.nu = np.array(wts)
        z = np.exp(2j*np.pi*self.t); self.w = z*(1-z)
        # full support for accurate rc
        self.tf = t[sup]; self.nuf = nu[sup]; zf = np.exp(2j*np.pi*self.tf); self.wf = zf*(1-zf)
        mp.mp.prec = 300
        self.wmp = []
        for tt in self.t:
            z = mp.expjpi(2*mp.mpf(tt)); self.wmp.append(z*(1-z))
        self.powcache = {}
    def pows(self, d):
        if d not in self.powcache:
            P = []
            for w in self.wmp:
                row = [mp.mpc(1)]
                for k in range(d): row.append(row[-1]*w)
                P.append(row)
            self.powcache[d] = P
        return self.powcache[d]
    def rc_fast(self, coeffs_desc):
        v = np.polyval(np.array(coeffs_desc, float), self.wf)
        return self.nuf @ np.log(np.abs(v))

def lattice(sp, d, weights, lam, bits=120):
    """weights: per support point (subsampled). Returns LLL-reduced basis rows (ascending coeff ints)."""
    P = sp.pows(d)
    S = mp.mpf(2)**bits
    cols = []
    for j in range(len(sp.t)):
        sc = mp.sqrt(mp.mpf(weights[j]))
        re = [int(mp.nint((P[j][k]*sc).real*S)) for k in range(d+1)]
        im = [int(mp.nint((P[j][k]*sc).imag*S)) for k in range(d+1)]
        cols.append(re); cols.append(im)
    lamI = max(1, int(mp.nint(mp.mpf(lam)*S)))
    M = [[lamI if i == k else 0 for i in range(d+1)] + [c[k] for c in cols] for k in range(d+1)]
    L = fmpz_mat(M).lll()
    out = []
    for r in range(d+1):
        v = [int(L[r, k]) for k in range(d+1)]
        if any(x % lamI for x in v): continue
        out.append([x//lamI for x in v])
    return out
