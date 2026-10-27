"""Fast numerics for C82 Doche construction: factor log-modulus on the circle, Mahler measures, active-set LP."""
import numpy as np, json
from fractions import Fraction
from flint import fmpz_poly, ctx
from scipy.optimize import linprog
import scipy.sparse as sp

def roots(coeffs_desc):
    ctx.prec = 200
    p = fmpz_poly([int(c) for c in coeffs_desc[::-1]])
    out = []
    for r, m in p.complex_roots():
        out += [complex(r)] * m
    return np.array(out)

def logabs(coeffs_desc, t, rts=None):
    if rts is None: rts = roots(coeffs_desc)
    w = np.exp(2j*np.pi*t); w = w*(1-w)
    out = np.full(len(t), np.log(abs(coeffs_desc[0])), float)
    for r in rts: out += np.log(np.abs(w - r))
    return out

def logmahlerT(coeffs_desc, rts=None):
    if rts is None: rts = roots(coeffs_desc)
    tot = np.log(abs(coeffs_desc[0]))
    for rho in rts:
        d = np.sqrt(1-4*rho+0j)
        for z in ((1+d)/2, (1-d)/2): tot += np.log(max(1.0, abs(z)))
    return tot

def grid(N):
    return (np.arange(N)+0.5)/N/2

class Lib:
    def __init__(self, den, num, t):
        self.t = t; self.den = den; self.num = num
        self.dR = np.array([len(r['coeffs'])-1 for r in den], float)
        self.dS = np.array([len(s['coeffs'])-1 for s in num], float)
        self.LM = np.array([logmahlerT(r['coeffs']) for r in den])
        self.LR = np.array([logabs(r['coeffs'], t) for r in den])
        self.LS = np.array([logabs(s['coeffs'], t) for s in num]) if num else np.zeros((0, len(t)))
    def add_num(self, s, ls=None):
        self.num.append(s); self.dS = np.append(self.dS, len(s['coeffs'])-1)
        if ls is None: ls = logabs(s['coeffs'], self.t)
        self.LS = np.vstack([self.LS, ls[None, :]])
    def add_den(self, r):
        self.den.append(r); self.dR = np.append(self.dR, len(r['coeffs'])-1)
        self.LM = np.append(self.LM, logmahlerT(r['coeffs']))
        self.LR = np.vstack([self.LR, logabs(r['coeffs'], self.t)[None, :]])
    def g(self, mu, w):
        return w @ self.LS - mu @ self.LR
    def value(self, mu, w):
        return self.LM @ mu + np.maximum(self.g(mu, w), 0).mean()

    def lp(self, mu0=None, w0=None, thr=0.02, budget=1-1e-7, verbose=False, maxit=30):
        """min LM.mu + mean(max(0,g)) s.t. dR.mu=1, dS.w<=budget. Active-set on grid points."""
        J, I, N = len(self.den), len(self.num), len(self.t)
        if mu0 is None:
            act = np.arange(N)
        else:
            gg = self.g(mu0, w0); act = np.flatnonzero(gg > -thr)
        for it in range(maxit):
            n = len(act)
            c = np.concatenate([self.LM, np.zeros(I), np.full(n, 1.0/N)])
            A1 = sp.hstack([sp.csr_matrix(-self.LR[:, act].T), sp.csr_matrix(self.LS[:, act].T), -sp.eye(n)])
            A2 = sp.csr_matrix(np.concatenate([np.zeros(J), self.dS, np.zeros(n)])[None, :])
            A_ub = sp.vstack([A1, A2]).tocsc()
            b_ub = np.concatenate([np.zeros(n), [budget]])
            A_eq = sp.csr_matrix(np.concatenate([self.dR, np.zeros(I+n)])[None, :])
            res = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=[1.0],
                          bounds=[(0, None)]*(J+I+n), method='highs-ds')
            assert res.status == 0, res.message
            mu, w = res.x[:J], res.x[J:J+I]
            gg = self.g(mu, w)
            viol = np.setdiff1d(np.flatnonzero(gg > -thr/4), act)
            if verbose: print(f'  it{it} n={n} lp={res.fun:.13f} true={self.value(mu,w):.13f} newpts={len(viol)}', flush=True)
            if len(viol) == 0: break
            act = np.union1d(act, np.flatnonzero(gg > -thr))
        dual = np.zeros(N); dual[act] = -res.ineqlin.marginals[:n]
        price = -res.ineqlin.marginals[n]; lam = res.eqlin.marginals[0]
        self.mu, self.w, self.dual, self.price, self.lam = mu, w, dual, price, lam
        return res.fun, mu, w

    def rc_num(self, ls, deg):
        return self.dual @ ls + self.price*deg
    def rc_den(self, lm, lr, deg):
        return lm - self.dual @ lr - self.lam*deg

def load(path):
    c = json.load(open(path)); return c, c['denominator'], c['numerator']
