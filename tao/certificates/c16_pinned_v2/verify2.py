"""Verify an upper-bound certificate for C16 (Brezis-Gallouet-Wainger remainder constant on T^2 = [0, 2pi]^2).

Claim: for every zero-mean u in H^2(T^2) (so delta(u) >= 1),
    4 pi ||u||_inf^2 / ||grad u||^2 - log delta(u) - log(1 + log delta(u)) <= T.

Proof structure (each step checked below in ball arithmetic):
 (1) Cauchy-Schwarz with weights |k|^2 (1 + mu |k|^2), mu > 0, and translation invariance:
        4 pi ||u||_inf^2 / ||grad u||^2 <= G(mu) (1 + mu delta) / pi,   G(mu) = sum_{k != 0} 1/(|k|^2 (1 + mu |k|^2)).
 (2) U(delta, mu) = G(mu)(1 + mu delta)/pi - log delta - log(1 + log delta) is convex in delta on [1, inf)
     (linear term; -log convex; log(1 + log delta) concave), so on [lo, hi] it is <= max(U(lo), U(hi)).
     The certificate's intervals cover [1, 2^18].
 (3) delta >= 2^18: take mu = 1/(lam delta), lam = log delta. For mu <= mu0 = 2^-20, lattice-point counting
     (#{k != 0 : |k| <= R} <= pi (R + 1/sqrt 2)^2) gives
        G(mu) <= G(mu0) + pi [log(mu0/mu) + pi c / a0 + c^2 / A0],   c = 1/sqrt 2, A0 = 1/mu0, a0 = sqrt A0,
     hence U <= 1 + B + (B + log lam)/lam - log(1 + 1/lam) <= 1 + B + (B + log Lam)/Lam for lam >= Lam = log 2^18,
     with B = G(mu0)/pi - log A0 + pi c/a0 + c^2/A0 (the last step needs (B + log lam)/lam decreasing: log Lam > 1 - B).
G(mu) is enclosed by exact row sums (pi/s) coth(pi s) and a trapezoid/midpoint bracket for rows x >= N (gsum.py).

    python verify.py certificate.json [--crosscheck K]
"""
import argparse, json, math
from flint import arb, ctx
from gsum import G_arb

ap = argparse.ArgumentParser(); ap.add_argument("cert"); ap.add_argument("--crosscheck", type=int, default=0)
ap.add_argument('--N', type=int, default=200000); ap.add_argument('--prec', type=int, default=192)
a = ap.parse_args(); ctx.prec = a.prec
c = json.load(open(a.cert)); T = arb(c["T"]); iv = [tuple(float.fromhex(x) for x in t) for t in c["intervals"]]
pi = arb.pi()

# coverage of [1, 2^18]
assert iv[0][0] == 1.0 and iv[-1][1] == 2.0 ** 18
for (lo, hi, mu), nxt in zip(iv, iv[1:] + [None]):
    assert lo < hi and mu > 0 and math.isfinite(mu)
    if nxt: assert hi == nxt[0], (hi, nxt[0])

Gc = {}; worst = None
for lo, hi, mu in iv:
    if mu not in Gc: Gc[mu] = G_arb(mu, N=a.N, prec=a.prec)
    for d in (lo, hi):
        d = arb(d); U = Gc[mu] / pi * (1 + arb(mu) * d) - d.log() - (1 + d.log()).log()
        assert U < T, (lo, hi, mu, U)
        worst = U if worst is None or U.upper() > worst.upper() else worst
print(f"[1, 2^18]: {len(iv)} intervals, {len(Gc)} values of mu; max endpoint bound {worst}  (N = {a.N}, {a.prec} bits)")

# tail delta >= 2^18
mu0 = 2.0 ** -20; A0 = arb(2) ** 20; a0 = A0.sqrt(); cc = 1 / arb(2).sqrt()
B = G_arb(mu0, N=a.N, prec=a.prec) / pi - A0.log() + pi * cc / a0 + cc * cc / A0
Lam = arb(2 ** 18).log()
assert Lam.log() > 1 - B
assert 1 / (Lam * 2 ** 18) <= arb(mu0)          # mu = 1/(lam delta) <= mu0 on the tail
tail = 1 + B + (B + Lam.log()) / Lam
assert tail < T
print(f"delta >= 2^18: bound {tail}  (B = {B})")
print(f"CERTIFIED: C16 <= {c['T']}")

if a.crosscheck:   # sanity check (not part of the proof): brute-force lattice sum over |k| <= R plus integral tail
    import numpy as np
    R = 3000; k = np.arange(-R, R + 1, dtype=float); r = k[:, None] ** 2 + k[None, :] ** 2; r = r[(r > 0) & (r <= R * R)]
    def G_brute(mu): return np.sum(1 / (r * (1 + mu * r))) + np.pi * np.log1p(1 / (mu * R * R))
    keys = sorted(Gc)[:: max(1, len(Gc) // a.crosscheck)]
    rel = max(abs(float(Gc[m].mid()) - G_brute(m)) / float(Gc[m].mid()) for m in keys)
    print(f"brute-force lattice-sum cross-check on {len(keys)} mu values: max relative deviation {rel:.1e}")
