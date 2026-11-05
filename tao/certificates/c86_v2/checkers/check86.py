"""Rigorous checker for C86 upper-bound certificates (piecewise-linear density + auxiliary polynomials).

Certificate (R. Pinho's format): exact rational edges e_0 < ... < e_n, node densities h_i >= 0 with h_0 = h_n = 0
(linear between nodes, zero outside), auxiliary polynomials Q_j (integer, constant first) with weights b_j >= 0,
and bound = T/m.  With mu the density, m = mu(R), T = int x dmu, U(z) = int log|z-x| dmu(x),
f(x) = U(x) - sum_j b_j log|Q_j(x)|, Smith's criterion gives C86 <= T/m when (Pinho, PROOF.md):
  (1) edges rational, ordered, positive; densities >= 0; Q_j primitive, irreducible, positive leading coefficient,
      pairwise distinct; b_j >= 0;
  (2) sum_j b_j deg Q_j < m;
  (3) int log|Q_j| dmu > 0 for every Q_j with b_j > 0;
  (4) f > 0 on every support cell, endpoints included;
  (5) the declared bound equals T/m exactly.
Potential: with slope jumps J_k of the density, U(z) = Re sum_k J_k G(z - e_k), G(w) = w^2 (log w / 2 - 3/4);
for nonreal z the principal log is continuous along the real segment, since Im(z - e_k) = Im z is fixed.
(4) uses, on [c - r, c + r], f(c + h) = f(c) + f'(c) h + R(h) with rigorous enclosures of R:
  knot with |t| > 2r (t = c - e_k):  G(t+h) - G(t) - h G'(t) = h^2/2 G''(xi),  G''(v) = log|v|, |v| in [|t|-r, |t|+r];
  knot with |t| <= 2r, 3r < 1:       G(t+h) - G(t) - h G'(t) in [r^2 (log(r/2) - 1), 0];
  root with d = |c - alpha| > 2r:    log|c+h-alpha| - log d - h Re 1/(c-alpha) in [-h^2/(2(d-r)^2), h^2/(2(d-r)^2)];
  root with d <= 2r:                 -log|x - alpha| >= -log(d + r)   (no linear term).
All arithmetic in Arb (python-flint) at 128 bits; decisions on ball bounds; intervals are bisected until accepted.

    python check86.py cert.json [--procs N] [--maxdepth D] [--beat X]
"""
import argparse, json, time
from fractions import Fraction as F
from multiprocessing import Pool
from flint import arb, acb, fmpz_poly, fmpq, ctx

PREC = 128
W = {}


def q2arb(x):
    return arb(fmpq(x.numerator, x.denominator))


def G(t):                      # real t (arb), G(0) = 0
    if t.contains(0):
        a = abs(t).upper(); a = arb(a)
        if a < arb(1) / 3:
            # |G(t)| <= a^2 (|log a| / 2 + 3/4) on [-a, a]
            bnd = a * a * (-a.log() / 2 + arb(3) / 4)
            return arb(0, bnd.upper())
        raise ValueError('wide ball at a knot')
    return t * t * (abs(t).log() / 2 - arb(3) / 4)


def Gp(t):                     # G'(t) = t (log|t| - 1), G'(0) = 0
    if t.contains(0):
        a = arb(abs(t).upper())
        return arb(0, (a * (1 - a.log())).upper())
    return t * (abs(t).log() - 1)


def Gc(w):                     # Re G(w) for complex w; on the real axis Re G is branch-independent
    a = abs(w)
    if a.contains(0) or a.lower() == 0:
        au = arb(a.upper())
        if au < arb(1) / 3:
            return arb(0, (au * au * (-au.log() / 2 + arb(3) / 4)).upper())
        raise ValueError('root ball too wide at a knot')
    return (w * w * (w.log() / 2 - acb(arb(3) / 4))).real


def init(cert):
    ctx.prec = PREC
    W['e'] = [q2arb(x) for x in cert['e']]; W['J'] = [q2arb(x) for x in cert['J']]
    W['eq'] = cert['e']
    W['roots'] = []; W['logc'] = arb(0)
    for Q, bj in cert['polys']:
        rr = fmpz_poly(Q).complex_roots()
        W['roots'].append([(z, bj) for z, _ in rr])
        W['logc'] += q2arb(bj) * arb(int(Q[-1])).log()


def cell_lower(c, r):
    """Lower bound of f on [c - r, c + r] (c, r exact rationals): returns arb."""
    ca, ra = q2arb(c), q2arb(r)
    f0 = arb(0); f1 = arb(0); rem_lo = arb(0); S2 = arb(0)
    r2 = ra * ra
    near_lo = r2 * ((ra / 2).log() - 1)
    for e, ek, J in zip(W['e'], W['eq'], W['J']):
        t = ca - e
        at = abs(c - ek)
        f0 += J * G(t); f1 += J * Gp(t)
        if at > 2 * r:
            # G''(xi) = log|xi| with |xi| in [|t|-r, |t|+r]; summed with signs, the remainder is h^2/2 * S2
            lo = (q2arb(at - r)).log(); hi = (q2arb(at + r)).log()
            S2 += J * arb((lo.lower() + hi.upper()) / 2, (hi.upper() - lo.lower()) / 2)
        else:
            if 3 * r >= 1: raise ValueError('interval too wide near a knot')
            if J > 0: rem_lo += J * near_lo
    rem_lo += arb(min(0, (S2 * r2 / 2).lower()))
    for rs in W['roots']:
        for z, bj in rs:
            b = q2arb(bj)
            w = acb(ca) - z
            d = abs(w)
            if d.lower() > 2 * ra:
                f0 -= b * d.log(); f1 -= b * (1 / w).real
                rem_lo -= b * r2 / (2 * (d - ra) ** 2)
            else:
                rem_lo -= b * (d + ra).log()
    f0 -= W['logc']
    return f0 - abs(f1) * ra + rem_lo


def check_cell(args):
    a, b, maxdepth = args
    stack = [(a, b, 0)]; nint = 0; worst = None
    while stack:
        lo, hi, dep = stack.pop(); nint += 1
        c, r = (lo + hi) / 2, (hi - lo) / 2
        try:
            L = cell_lower(c, r)
            ok = L > 0
        except ValueError:
            ok = False; L = None
        if ok:
            if worst is None or L.lower() < worst: worst = float(L.lower())
            continue
        if dep >= maxdepth:
            return (str(a), str(b), False, nint, None)
        stack += [(lo, c, dep + 1), (c, hi, dep + 1)]
    return (str(a), str(b), True, nint, worst)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('cert'); ap.add_argument('--procs', type=int, default=8)
    ap.add_argument('--maxdepth', type=int, default=22); ap.add_argument('--beat', default=None)
    ap.add_argument('--out', default=None); ap.add_argument('--limit', type=int, default=0)
    a = ap.parse_args(); ctx.prec = PREC; t0 = time.time()
    c = json.load(open(a.cert))
    e = [F(x) for x in c['edges']]; h = [F(x) for x in c['node_densities']]
    A = c['A']; B = [F(x) for x in c['b']]
    n = len(e)
    # (1) structure
    assert all(e[i] < e[i + 1] for i in range(n - 1)) and e[0] > 0
    assert all(x >= 0 for x in h) and h[0] == 0 and h[-1] == 0
    assert all(bj >= 0 for bj in B) and len(A) == len(B)
    seen = set()
    for Q in A:
        P = fmpz_poly(Q)
        assert P.degree() >= 1 and int(P.coeffs()[-1]) > 0
        from math import gcd
        g = 0
        for t in Q: g = gcd(g, abs(int(t)))
        assert g == 1, 'not primitive'
        fac = P.factor()[1]
        assert len(fac) == 1 and fac[0][1] == 1, 'not irreducible'
        assert tuple(Q) not in seen; seen.add(tuple(Q))
    # (5) exact mass, moment, bound
    m = sum((e[i + 1] - e[i]) * (h[i] + h[i + 1]) / 2 for i in range(n - 1))
    T = sum((e[i + 1] - e[i]) * (h[i] * (2 * e[i] + e[i + 1]) + h[i + 1] * (e[i] + 2 * e[i + 1])) / 6 for i in range(n - 1))
    bound = T / m
    assert bound == F(c['bound']), 'declared bound differs from T/m'
    # (2) degree condition
    assert sum(bj * (len(Q) - 1) for Q, bj in zip(A, B)) < m
    # slope jumps and identities
    s = [F(0)] + [(h[i + 1] - h[i]) / (e[i + 1] - e[i]) for i in range(n - 1)] + [F(0)]
    J = [s[i + 1] - s[i] for i in range(n)]
    assert sum(J) == 0 and sum(j * x for j, x in zip(J, e)) == 0 and sum(j * x * x for j, x in zip(J, e)) / 2 == m
    keepk = [i for i in range(n) if J[i] != 0]
    ek = [e[i] for i in keepk]; Jk = [J[i] for i in keepk]
    ea = [q2arb(x) for x in ek]; Ja = [q2arb(x) for x in Jk]
    # roots of the selected polynomials (certified balls), simple
    roots = []
    logc = arb(0)
    for Q, bj in zip(A, B):
        if bj == 0: continue
        rr = fmpz_poly(Q).complex_roots()
        assert all(mult == 1 for _, mult in rr) and len(rr) == len(Q) - 1
        roots.append(([z for z, _ in rr], bj))
        logc += q2arb(bj) * arb(int(Q[-1])).log()
    # (3) integrals int log|Q| dmu = m log lc + sum_roots U(alpha)
    mins = None
    for (rs, bj), Q in zip(roots, [Q for Q, bj in zip(A, B) if bj > 0]):
        I = q2arb(m) * arb(int(Q[-1])).log()
        for z in rs:
            I += sum((Jv * Gc(z - acb(ev)) for Jv, ev in zip(Ja, ea)), arb(0))
        assert I > 0, f'integral not certified positive for {Q[:4]}...: {I}'
        mins = I if mins is None or I.lower() < mins.lower() else mins
    print(f'structure, bound {float(bound):.13f}, degree condition: OK; min integral {mins}  ({time.time() - t0:.0f}s)', flush=True)
    # (4) positivity on support cells
    sup = [i for i in range(n - 1) if h[i] > 0 or h[i + 1] > 0]
    if a.limit: sup = sup[::max(1, len(sup) // a.limit)][:a.limit]
    cert = {'e': ek, 'J': Jk, 'polys': [(Q, bj) for Q, bj in zip(A, B) if bj > 0]}
    with Pool(a.procs, initializer=init, initargs=(cert,)) as pool:
        res = pool.map(check_cell, [(e[i], e[i + 1], a.maxdepth) for i in sup], chunksize=8)
    bad = [x for x in res if not x[2]]
    worst = min(x[4] for x in res if x[2])
    nint = sum(x[3] for x in res)
    print(f'positivity: {len(sup)} support cells, {nint} intervals, {len(bad)} failed; min certified lower bound {worst:.3e}  ({time.time() - t0:.0f}s)', flush=True)
    ok = not bad
    if a.beat: print(f'bound {float(bound):.13f} < {a.beat}: {bound < F(a.beat)}')
    print('VALID' if ok else 'NOT CERTIFIED')
    if a.out:
        json.dump({'valid': ok, 'bound': str(bound), 'bound_float': float(bound), 'support_cells': len(sup), 'intervals': nint,
                   'failed': bad[:20], 'min_positivity_lower_bound': worst, 'seconds': time.time() - t0}, open(a.out, 'w'), indent=1)


if __name__ == '__main__':
    main()
