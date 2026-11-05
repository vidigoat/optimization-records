"""Fine upper-bound certificate for C16 at target T (see build.py / verify.py for the inequality).
mu per interval [lo, hi]: the balance point mu G(mu) = pi (Lam(hi) - Lam(lo)) / (hi - lo), Lam = log d + log(1 + log d),
clipped to [mu*(lo), mu*(hi)] where mu*(d) solves G'(mu)(1 + mu d) + G(mu) d = 0; G and G' = -A^2 H in double
(only the choice of mu, which is then checked rigorously).  Rigorous check: Arb, N = 200000 rows, 192 bits.
    python build2.py T out.json [procs]"""
import json, math, sys
import numpy as np
from multiprocessing import Pool
from gsum import G_arb
from flint import arb, ctx

NR = 200000
x = np.arange(1, 400000, dtype=float)
def GH(mu):
    A = 1.0 / mu; a = math.sqrt(A); s = np.sqrt(x * x + A); pi = math.pi
    rows = pi / x / np.tanh(pi * x) - pi / s / np.tanh(pi * s)
    M = x[-1] + 0.5
    G = pi ** 2 / 3 - (pi / a / math.tanh(pi * a) - 1 / A) + 2 * (np.sum(rows[::-1]) + pi * math.log((M + math.sqrt(M * M + A)) / (2 * M)))
    def rs(s): return pi / (2 * s ** 3) / np.tanh(pi * s) + pi * pi / (2 * s * s) / np.sinh(np.minimum(pi * s, 700)) ** 2
    H = rs(a) - 1 / (A * A) + 2 * (np.sum(rs(s)[::-1]) + pi / (2 * A) * (1 - M / math.sqrt(M * M + A)))
    return G, -A * A * H
def bis(f, lo, hi, it=200):
    flo = f(lo)
    for _ in range(it):
        m = math.sqrt(lo * hi); fm = f(m)
        if (fm > 0) == (flo > 0): lo, flo = m, fm
        else: hi = m
        if hi / lo - 1 < 1e-15: break
    return math.sqrt(lo * hi)
def mustar(d):
    def F(mu): G, Gp = GH(mu); return Gp * (1 + mu * d) + G * d
    return bis(F, 1e-9, 1e3)
def Lam(d): return math.log(d) + math.log1p(math.log(d))
def choose_mu(lo, hi):
    m1, m2 = mustar(lo), mustar(hi)
    import mpmath as mp
    mp.mp.dps = 50
    L = lambda d: mp.log(d) + mp.log(1 + mp.log(d))
    t = float(mp.pi * (L(mp.mpf(hi)) - L(mp.mpf(lo))) / (mp.mpf(hi) - mp.mpf(lo)))
    mb = bis(lambda mu: mu * GH(mu)[0] - t, 1e-9, 1e3)
    return min(max(mb, min(m1, m2)), max(m1, m2))
TS = None
def init(s):
    global TS
    TS = s
def work(iv):
    lo, hi = iv; ctx.prec = 192
    mu = choose_mu(lo, hi); G = G_arb(mu, N=NR, prec=192); pi = arb.pi()
    us = []
    for d in (lo, hi):
        d = arb(d); us.append(G / pi * (1 + arb(mu) * d) - d.log() - (1 + d.log()).log())
    ok = all(u < arb(TS) for u in us)
    return lo, hi, mu, ok, max(float(u.upper()) for u in us), str(max(us, key=lambda u: u.upper()))
if __name__ == '__main__':
    TS = sys.argv[1]; out = sys.argv[2]; procs = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    c = json.load(open('certificate.json'))
    todo = [(float.fromhex(a), float.fromhex(b)) for a, b, _ in c['intervals']]
    done = []; worst = None; rnd = 0
    with Pool(procs, initializer=init, initargs=(TS,)) as pool:
        while todo:
            res = pool.map(work, todo, chunksize=4); todo = []
            for lo, hi, mu, ok, ub, us in res:
                if ok:
                    done.append((lo.hex(), hi.hex(), mu.hex()))
                    if worst is None or ub > worst[0]: worst = (ub, us)
                else:
                    m = (lo + hi) / 2
                    if not (lo < m < hi): raise SystemExit(f'cannot split {lo} {hi}')
                    todo += [(lo, m), (m, hi)]
            rnd += 1; print(f'round {rnd}: {len(done)} accepted, {len(todo)} to split', flush=True)
    done.sort(key=lambda t: float.fromhex(t[0]))
    json.dump({'T': TS, 'intervals': done}, open(out, 'w'))
    print(f'{len(done)} intervals; worst endpoint bound {worst[1]}')
