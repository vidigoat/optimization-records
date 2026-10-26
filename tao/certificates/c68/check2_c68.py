#!/usr/bin/env python3
"""Second, separately written checker (mpmath interval arithmetic, no python-flint).

Same mathematical hypotheses as check_c68.py, different arithmetic library and
code.  Decisions use interval endpoints only: a quantity is accepted as positive
iff the LOWER endpoint of its enclosing interval is > 0.
"""
import json, sys, argparse, time
from fractions import Fraction as Fr
from mpmath import iv, mpf

if not __debug__:
    raise SystemExit('refusing to run with -O')


def I(x):
    x = Fr(x)
    return iv.mpf(x.numerator) / iv.mpf(x.denominator)


def lo(x):
    return x.a


def hi(x):
    return x.b


def pos(x):
    return lo(x) > 0


def sinh(x):
    return (iv.exp(x) - iv.exp(-x)) / 2


def cosh(x):
    return (iv.exp(x) + iv.exp(-x)) / 2


def main(path, prec):
    iv.prec = prec
    d = json.load(open(path))
    model = d.get('model') or d.get('analytic_model')
    c0, R0 = Fr(d['c']), Fr(d['R'])
    ts = [Fr(x) for x in d['nodes']]; ws = [Fr(x) for x in d['weights']]
    # exact structural hypotheses
    assert 0 < c0 < 1 and c0 * c0 < R0 < 1 and R0 >= c0
    assert ts[0] == c0 * c0 and ts[-1] == R0
    assert all(a < b for a, b in zip(ts, ts[1:]))
    assert all(x > 0 for x in ws) and sum(ws) == R0
    pi = iv.pi
    w = -iv.log(I(c0))
    R = I(R0)
    out = {'model': model, 'c': str(c0)}

    if model in ('two_covering_zeros_hyperbolic_jung', 'pinho_jung'):
        X = sinh(pi * pi / (2 * w)) ** 2

        def gl(t):   # lower bound of g = 1/G at atom t
            y = iv.sin(pi * iv.log(I(t)) / (2 * w)) ** 2 / X
            H = 3 * y * (2 + y)
            v = 2 * iv.sqrt(H) / (iv.sqrt(H) + iv.sqrt(H + 4))
            return lo(v)
        ell = Fr(d['linear_slope'])
        U = -pi * iv.log(R) / (2 * w)
        assert lo(U) > 0 and hi(U * U) < 6
        a = (pi / (2 * w)) ** 2 * (1 - U * U / 6) ** 2 / X
        z = a * (1 - R) ** 2
        bound = 2 * iv.sqrt(6 * a) / (iv.sqrt(3 * z * (2 + z)) + iv.sqrt(4 + 6 * z + 3 * z * z))
        assert lo(bound) > I(ell).b, 'slope'
        pieces = [(R0, Fr(1), I(ell))]
        sec = [(Fr(d['second_tail_t0']), Fr(d['second_tail_t1']))]
        Kfix = int(d['second_tail_k'])
    elif model == 'flux_halfdiameter':
        A = pi * pi / w
        q = iv.exp(-A)
        Psi = (1 + q) ** 2 / (8 * q)
        assert hi(1 / cosh(A / 2) ** 2) < mpf(1) / 3
        assert hi(8 * cosh(A / 6) ** 4) < lo(Psi)

        def gl(t):
            s = iv.sin(pi * iv.log(I(t)) / (2 * w)) ** 2
            v = 2 * iv.sqrt(2 * s) / (iv.sqrt(2 * s) + iv.sqrt(2 * s + Psi))
            return lo(v)
        tp = [Fr(x) for x in d['tail_partition']]
        assert tp[0] == R0 and tp[-1] == 1 and all(a < b for a, b in zip(tp, tp[1:]))
        pieces = []
        for a0, b0 in zip(tp, tp[1:]):
            dh, dl = 1 - I(a0), 1 - I(b0)
            U = -pi * iv.log(I(a0)) / (2 * w)
            assert lo(U) > 0 and hi(U * U) < 6
            bb = (pi / (2 * w)) * (1 + dl / 2) * (1 - U * U / 6)
            r2 = iv.sqrt(2)
            l = 2 * r2 * bb / (r2 * bb * dh + iv.sqrt(2 * bb * bb * dh * dh + Psi))
            assert lo(l) > 0
            pieces.append((a0, b0, iv.mpf(lo(l))))
        sec = [(a0, b0) for (a0, b0, _) in pieces if b0 < 1]
        Kfix = None
    else:
        raise SystemExit('unknown model')

    # beta upper bounds at atoms (beta = 1 - g, beta(c^2)=1)
    beta = []
    for t in ts:
        if t == c0 * c0:
            beta.append(iv.mpf(1))
        else:
            gv = gl(t)
            assert gv > 0
            beta.append(iv.mpf(1) - iv.mpf(gv))     # interval whose upper end >= true beta
    # first family tail index
    K1 = int(-(-R0 // ws[-1])) - 1                 # smallest K with w_last (K+1) >= R
    while ws[-1] * K1 >= R0: K1 -= 1
    while ws[-1] * (K1 + 1) < R0: K1 += 1
    assert ws[-1] * (K1 + 1) >= R0
    T = [I(t) for t in ts]; W = [I(x) for x in ws]
    BW = [wj * bj for wj, bj in zip(W, beta)]
    nodes_p = sorted(set([p[0] for p in pieces] + [p[1] for p in pieces]))
    PI = {p: I(p) for p in nodes_p}
    pp = dict(PI)                                  # u^{k+1}
    pw = [iv.mpf(1) for _ in T]
    Rp = R
    K2 = None; mA = mB = None
    k = 0
    while True:
        k += 1
        pw = [x * t for x, t in zip(pw, T)]
        Rp = Rp * R
        for p in nodes_p: pp[p] = pp[p] * PI[p]
        s = Rp / (k + 1)
        if k < K1:
            Ak = sum([wj * x for wj, x in zip(W, pw)], iv.mpf(0)) - s
            assert pos(Ak), ('A', k)
            r = lo(Ak / s); mA = r if mA is None or r < mA else mA
        if K2 is None:
            tot = iv.mpf(0)
            for (a0, b0, l) in pieces:
                tot += l * ((pp[b0] - pp[a0]) / (k + 1) - (pp[b0] * PI[b0] - pp[a0] * PI[a0]) / (k + 2))
            Bk = tot + s - sum([bw * x for bw, x in zip(BW, pw)], iv.mpf(0))
            assert pos(Bk), ('B', k)
            r = lo(Bk / s); mB = r if mB is None or r < mB else mB
            nk = k + 1
            if (Kfix is not None and nk == Kfix) or (Kfix is None and nk % 50 == 0):
                neg = sum([wj * x * t for wj, x, t in zip(W, pw, T)], iv.mpf(0))
                for (a0, b0) in sec:
                    lval = [p[2] for p in pieces if p[0] <= a0 and b0 <= p[1]][0]
                    posv = lval * (1 - I(b0)) * (I(b0) - I(a0)) * I(a0) ** nk
                    if lo(posv) > hi(neg) and a0 >= R0:
                        K2 = nk; out['B_tail'] = [nk, str(a0), str(b0)]
                        break
                if Kfix is not None:
                    assert K2 is not None
        if k >= K1 - 1 and K2 is not None:
            break
    out.update(status='PASS', A_finite_through=K1 - 1, B_finite_through=K2 - 1,
               minA=float(mA), minB=float(mB))
    return out


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('certificate'); ap.add_argument('--prec', type=int, default=256)
    a = ap.parse_args()
    t0 = time.time()
    r = main(a.certificate, a.prec)
    r['seconds'] = round(time.time() - t0, 1)
    print(json.dumps(r, indent=1))
