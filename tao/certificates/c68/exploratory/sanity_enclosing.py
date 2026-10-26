"""Numerical sanity tests (NOT part of the proof) of the analytic lemmas:
For annular Blaschke products omega (period condition satisfied), on circles |z|=r:
 (i)  pairwise pseudohyperbolic distances <= F_P(theta) * exp(-E),  E = sigma sech^2(A/2)
 (ii) minimal enclosing hyperbolic radius rho satisfies (1+cosh rho)/2 <= G_H(t)
"""
import numpy as np, mpmath as mp
from scipy.optimize import minimize
rng = np.random.default_rng(1)
mp.mp.dps = 30

def setup(c):
    w = -np.log(c); A = np.pi**2 / w; q = np.exp(-2 * A)
    return w, A, q

def Ecov(zeta, w):  # strip -> UHP
    return np.exp(1j * np.pi * (zeta + w) / w)

def blaschke_orbit(Z, a0, q, K=12):
    """prod over k of normalised UHP Blaschke factors with zeros a0 q^k."""
    out = np.ones_like(Z, dtype=complex)
    for k in range(-K, K + 1):
        al = a0 * q**k
        if k >= 0:
            out *= (Z - al) / (Z - np.conj(al))
        else:
            out *= (1 - Z / al) / (1 - Z / np.conj(al))
    return out

def make_omega(c, zeros):
    """zeros: list of complex points in the annulus. returns omega on annulus (callable)."""
    w, A, q = setup(c)
    a0s = [Ecov(np.log(z), w) for z in zeros]
    def om(z):
        Z = Ecov(np.log(z), w)
        v = np.ones_like(Z, dtype=complex)
        for a0 in a0s:
            v *= blaschke_orbit(Z, a0, q)
        return v
    return om

def pdist(a, b):
    return np.abs(a - b) / np.abs(1 - np.conj(b) * a)

def hdist(a, b):
    return 2 * np.arctanh(np.minimum(pdist(a, b), 1 - 1e-16))

def min_enclosing(P):
    """minimal enclosing hyperbolic radius of points P (Poincare disk)."""
    best = np.inf
    starts = [P[0], P[len(P) // 2], np.mean(P)]
    for s0 in starts:
        # parametrise centre via Mobius recentering for conditioning
        def f(x):
            cen = s0 + (x[0] + 1j * x[1]) * (1 - abs(s0)**2)
            if abs(cen) >= 1: return 1e9
            return np.max(hdist(P, cen))
        r = minimize(f, [0.0, 0.0], method='Nelder-Mead', options=dict(xatol=1e-12, fatol=1e-12, maxiter=20000))
        best = min(best, r.fun)
    return best

def bounds(c, r, theta):
    w, A, q = setup(c)
    t = r * r
    sig = np.sin(np.pi * np.log(t) / (2 * w))**2
    v = A * theta / np.pi
    S1 = np.sinh(v / 2)**2; S2 = np.sinh(A - v / 2)**2
    FP = np.sqrt(S1 / (S1 + sig) * S2 / (S2 + sig))
    E = sig / np.cosh(A / 2)**2
    qq = np.exp(-A); Psi = (1 + qq)**2 / (8 * qq)
    GH = (1 + np.sqrt(1 + Psi / (2 * sig))) / 2
    return FP * np.exp(-E), GH

if __name__ == "__main__":
  c = 0.446
  w, A, q = setup(c)
  worst_pair = 0; worst_G = 0
  for trial in range(40):
      r = rng.uniform(c * 1.02, 0.98)
      deg = rng.integers(2, 4)
      # choose zeros with sum of phi == integer (period condition)
      phis = rng.uniform(0.02, 0.98, deg - 1)
      last = (np.ceil(phis.sum()) - phis.sum())
      if last < 0.02: last += 1
      if last > 0.98: continue
      phis = np.r_[phis, last]
      if trial % 2 == 0:
          # include the extremal-type configuration: zero at r, partner at c/r antipodal-ish
          phis = np.array([1 + np.log(r) / w, -np.log(r) / w])
          args = np.array([0.0, rng.normal(0, 0.3)])
      else:
          args = rng.uniform(0, 2 * np.pi, len(phis))
      radii = np.exp(w * (phis - 1))
      zeros = radii * np.exp(1j * args)
      om = make_omega(c, zeros)
      # periodicity check
      z0 = r * np.exp(0.3j)
      Z = Ecov(np.log(z0), w)
      per = abs(om(np.exp(np.log(z0) + 2j * np.pi)) / om(z0) - 1) if False else 0
      th = np.linspace(0, 2 * np.pi, 721)[:-1]
      zz = r * np.exp(1j * th)
      vals = om(zz)
      # check periodicity numerically via lifted value at theta=2pi
      Zp = Ecov(np.log(r) + 2j * np.pi, w)
      v2 = np.ones(1, complex)
      for zz0 in zeros:
          v2 *= blaschke_orbit(np.array([Zp]), Ecov(np.log(zz0), w), q)
      perr = abs(v2[0] - vals[0])
      # (i) pairwise
      i0 = np.arange(0, 720, 7)
      ratio = 0
      for i in i0:
          dth = np.abs(th - th[i])
          b, GH = bounds(c, r, dth)
          d = pdist(vals, vals[i])
          mask = (b < 1 - 1e-14) & (dth > 1e-9)
          ratio = max(ratio, np.max(d[mask] / b[mask]))  # >1 means violation
      worst_pair = max(worst_pair, ratio)
      # (ii) enclosing
      rho = min_enclosing(vals[::4])
      _, GH = bounds(c, r, np.pi)
      Gt = (1 + np.cosh(rho)) / 2
      worst_G = max(worst_G, Gt / GH)
      print('trial %2d r=%.3f deg=%d period_err=%.1e  max dist/bound=%.12f  G_true/G_H=%.6f' % (trial, r, len(phis), perr, ratio, Gt / GH))
  print('WORST pair ratio (must be <=1):', worst_pair, ' WORST G ratio (must be <=1):', worst_G)
