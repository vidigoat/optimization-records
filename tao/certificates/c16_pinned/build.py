"""Build an upper-bound certificate C16 <= T: intervals [lo, hi] of delta covering [1, 2^18], each with a mu > 0.
For mu > 0 and any zero-mean u in H^2(T^2):  4 pi |u(0)|^2 / ||grad u||^2 <= G(mu) (1 + mu delta(u)) / pi  (Cauchy-Schwarz).
U(delta, mu) = G(mu)(1+mu delta)/pi - log delta - log(1 + log delta) is convex in delta >= 1, so on [lo, hi] it is
at most max(U(lo), U(hi)).  delta >= 2^18 is covered analytically in verify.py.
    python build.py T out.json"""
import json, math, sys
from scipy.optimize import minimize_scalar
from gsum import G_float, G_arb
from flint import arb, ctx

T = float(sys.argv[1]); out = sys.argv[2]; ctx.prec = 160
Gc = {}
def Gf(mu): return G_float(mu, N=20000)
def Uf(d, mu): return Gf(mu) / math.pi * (1 + mu * d) - math.log(d) - math.log1p(math.log(d))
def best_mu(lo, hi):
    f = lambda lm: max(Uf(lo, math.exp(lm)), Uf(hi, math.exp(lm)))
    r = minimize_scalar(f, bounds=(math.log(1e-8), math.log(1e3)), method="bounded", options={"xatol": 1e-10})
    return math.exp(r.x)
def Ua(d, mu):
    if mu not in Gc: Gc[mu] = G_arb(mu)
    d = arb(d); return Gc[mu] / arb.pi() * (1 + arb(mu) * d) - d.log() - (1 + d.log()).log()

cert = []; worst = -1e9
edges = [1.0] + [(2.0 ** 18) ** (i / 256) for i in range(1, 256)] + [2.0 ** 18]
stack = list(zip(edges[:-1], edges[1:]))
while stack:
    lo, hi = stack.pop()
    mu = best_mu(lo, hi)
    if max(Uf(lo, mu), Uf(hi, mu)) < T:
        b = max(Ua(lo, mu).upper(), Ua(hi, mu).upper())
        if b < T:
            cert.append((lo.hex(), hi.hex(), mu.hex())); worst = max(worst, float(b)); continue
    m = (lo + hi) / 2
    if m <= lo or m >= hi: raise SystemExit(f"cannot split [{lo}, {hi}]")
    stack += [(m, hi), (lo, m)]
cert.sort(key=lambda t: float.fromhex(t[0]))
json.dump({"T": repr(T), "intervals": cert}, open(out, "w"))
print(f"{len(cert)} intervals, {len(Gc)} mu values, worst endpoint upper bound {worst!r}")
