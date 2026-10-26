#!/usr/bin/env python3
"""Rigorous checker for moment certificates for Korenblum's constant (C68).

Independent of any optimiser.  All certificate data are exact rationals; every
inequality is decided by Arb ball arithmetic (python-flint) with outward
rounding.  Two analytic models are supported:

  pinho_jung        : Pinho (2026) two-covering-zero bound + hyperbolic Jung,
                      linear tail envelope l*(1-t) on [R,1]   (replay of 0.432)
  flux_halfdiameter : two-covering-zero bound sharpened by the flux/period term
                      E = sigma*sech^2(A/2), with the half-diameter enclosing
                      lemma:  G(t) <= (1+sqrt(1+Psi/(2 sigma(t))))/2,
                      Psi = (1+q)^2/(8q), q = exp(-A), A = pi^2/w, w = -log c.

Wikstrom's moment criterion: a probability measure mu on [c^2,1), no atom at 1,
with  int t^k dmu >= 1/(k+1) (k>=0)  and  int beta t^k dmu <= 1/(k+1) (k>=1),
beta = 1 - 1/G.  For mu = sum_j w_j delta_{t_j} + 1_[R,1] dt these become
   A_k = sum w_j t_j^k - R^{k+1}/(k+1) >= 0,
   B_k = int_R^1 g t^k dt + R^{k+1}/(k+1) - sum w_j beta_j t_j^k >= 0,  g = 1/G.
Every k is covered: finite range by direct enclosure, the rest by tail lemmas.
"""
import sys, json, argparse, hashlib, time
from fractions import Fraction as Fr
from flint import arb, fmpq, ctx

if not __debug__:
    raise SystemExit('refusing to run with -O: assertions are the proof')


def A_(x):
    """Exact rational -> arb ball (tight enclosure)."""
    x = Fr(x)
    return arb(fmpq(x.numerator, x.denominator))


def is_pos(x):
    return x > 0          # arb comparison: True only if the whole ball is > 0


def common_checks(d):
    c, R = Fr(d['c']), Fr(d['R'])
    ts = [Fr(x) for x in d['nodes']]
    ws = [Fr(x) for x in d['weights']]
    assert 0 < c < 1 and c * c < R < 1
    assert len(ts) == len(ws) >= 2
    assert ts[0] == c * c, 'first node must be c^2'
    assert ts[-1] == R, 'last node must be R'
    assert all(x < y for x, y in zip(ts, ts[1:])), 'nodes strictly increasing'
    assert all(v > 0 for v in ws), 'weights positive'
    assert sum(ws) == R, 'atomic mass must equal R exactly (probability measure)'
    assert R >= c, 'need R >= c for monotonicity of sigma on [R,1]'
    return c, R, ts, ws


def run(d, prec_bits=320, verbose=True):
    ctx.prec = prec_bits
    t0 = time.time()
    model = d.get('model') or d.get('analytic_model')
    c0, R0, ts, ws = common_checks(d)
    c, R = A_(c0), A_(R0)
    pi = arb.pi()
    w = -c.log()
    T = [A_(x) for x in ts]
    W = [A_(x) for x in ws]
    n = len(T)
    info = {'model': model, 'c': str(c0), 'R': str(R0), 'atoms': n}

    # ---------------- analytic model: beta at atoms, tail pieces -------------
    if model == 'two_covering_zeros_hyperbolic_jung' or model == 'pinho_jung':
        X = (pi * pi / (2 * w)).sinh() ** 2

        def g_atom(t):
            y = (pi * t.log() / (2 * w)).sin() ** 2 / X
            H = 3 * y * (2 + y)
            return 2 * H.sqrt() / (H.sqrt() + (H + 4).sqrt())
        l0 = Fr(d['linear_slope'])
        U = -pi * R.log() / (2 * w)
        assert U > 0 and U * U < 6
        aa = (pi / (2 * w)) ** 2 * (1 - U * U / 6) ** 2 / X
        z = aa * (1 - R) ** 2
        slope = 2 * (6 * aa).sqrt() / ((3 * z * (2 + z)).sqrt() + (4 + 6 * z + 3 * z * z).sqrt())
        assert slope > A_(l0), 'linear envelope not below derived bound'
        pieces = [(R0, Fr(1), A_(l0))]           # g >= l0 (1-t) on [R,1]
        info['linear_slope_bound_lower'] = slope.lower().str(20)
        tail_pieces_for_K2 = [(Fr(d['second_tail_t0']), Fr(d['second_tail_t1']), A_(l0))]
        assert R0 <= tail_pieces_for_K2[0][0] < tail_pieces_for_K2[0][1] < 1
        K2_given = int(d['second_tail_k'])
    elif model == 'flux_halfdiameter':
        A = pi * pi / w
        q = (-A).exp()
        Psi = (1 + q) ** 2 / (8 * q)
        # hypotheses of the analytic lemmas
        sech2 = 1 / (A / 2).cosh() ** 2
        assert sech2 < arb(1) / 3, 'convexity hypothesis sech^2(A/2) <= 1/3'
        assert 8 * (A / 6).cosh() ** 4 < Psi, 'simplex hypothesis 8 cosh^4(A/6) <= Psi'
        info['A'] = A.str(15); info['Psi'] = Psi.str(15)

        def g_atom(t):
            s = (pi * t.log() / (2 * w)).sin() ** 2
            return 2 * (2 * s).sqrt() / ((2 * s).sqrt() + (2 * s + Psi).sqrt())
        tp = [Fr(x) for x in d['tail_partition']]
        assert tp[0] == R0 and tp[-1] == 1 and all(x < y for x, y in zip(tp, tp[1:]))
        pieces = []
        for a0, b0 in zip(tp, tp[1:]):
            a, b = A_(a0), A_(b0)
            dhi, dlo = 1 - a, 1 - b
            U = -pi * a.log() / (2 * w)
            assert U > 0 and U * U < 6
            bb = (pi / (2 * w)) * (1 + dlo / 2) * (1 - U * U / 6)
            assert bb > 0
            l = 2 * arb(2).sqrt() * bb / (arb(2).sqrt() * bb * dhi + (2 * bb * bb * dhi * dhi + Psi).sqrt())
            assert l > 0
            pieces.append((a0, b0, arb(l.lower())))     # g >= l (1-t) on [a,b]
        tail_pieces_for_K2 = [p for p in pieces if p[1] < 1]
        K2_given = None
    else:
        raise SystemExit('unknown model %r' % model)

    beta = []
    for t0_, t in zip(ts, T):
        if t0_ == c0 * c0:
            beta.append(arb(1))                      # convention beta(c^2)=1
        else:
            g = g_atom(t)
            assert g > 0 and g < 1
            b = 1 - g
            beta.append(arb(b.upper()) if b.upper() < 1 else arb(1))
    info['beta_upper'] = [b.str(12) for b in beta]

    # ---------------- first tail: A_k for k >= K1 by the atom at R ----------
    K1 = 1
    while ws[-1] * (K1 + 1) < R0:
        K1 += 1
    assert ws[-1] * (K1 + 1) >= R0
    info['A_tail_from'] = K1

    # ---------------- finite checks ------------------------------------------
    pw = [arb(1)] * n                       # t_j^k
    pts = sorted(set([p[0] for p in pieces] + [p[1] for p in pieces]))
    PA = {p: A_(p) for p in pts}
    ppow = {p: PA[p] for p in pts}          # u^{k+1}
    Rk1 = R                                 # R^{k+1}
    minA = minB = None; kA = kB = 0
    K2 = None
    kmax_guard = 10**6
    k = 0
    while True:
        k += 1
        assert k < kmax_guard
        pw = [x * t for x, t in zip(pw, T)]
        Rk1 = Rk1 * R
        for p in pts:
            ppow[p] = ppow[p] * PA[p]
        sk = Rk1 / (k + 1)
        if k < K1:
            Ak = sum(wj * x for wj, x in zip(W, pw)) - sk
            assert is_pos(Ak), ('A_k fails', k, Ak.str(10))
            rel = (Ak / sk).lower()
            if minA is None or rel < minA: minA, kA = rel, k
        if K2 is None:
            # J(u) = u^{k+1}/(k+1) - u^{k+2}/(k+2)
            I = arb(0)
            for (a0, b0, l) in pieces:
                Jb = ppow[b0] / (k + 1) - ppow[b0] * PA[b0] / (k + 2)
                Ja = ppow[a0] / (k + 1) - ppow[a0] * PA[a0] / (k + 2)
                I += l * (Jb - Ja)
            Bk = I + sk - sum(wj * bj * x for wj, bj, x in zip(W, beta, pw))
            assert is_pos(Bk), ('B_k fails', k, Bk.str(10))
            rel = (Bk / sk).lower()
            if minB is None or rel < minB: minB, kB = rel, k
            # try to close the second family by the tail lemma
            want = (K2_given is not None and k + 1 == K2_given) or (K2_given is None and (k + 1) % 50 == 0)
            if want:
                kk = k + 1
                neg = sum(wj * x * t for wj, x, t in zip(W, pw, T))       # sum w t^{kk}
                for (a0, b0, l) in tail_pieces_for_K2:
                    a, b = A_(a0), A_(b0)
                    pos = l * (1 - b) * (b - a) * a ** kk
                    if pos > neg:
                        assert a0 >= R0
                        K2 = kk
                        info['B_tail_from'] = kk
                        info['B_tail_piece'] = [str(a0), str(b0)]
                        info['B_tail_gap_lower'] = (pos - neg).lower().str(10)
                        break
                if K2_given is not None:
                    assert K2 is not None, 'given second tail fails'
        if k >= K1 - 1 and K2 is not None:
            break
    info.update({'status': 'PASS', 'A_finite_through': K1 - 1, 'B_finite_through': K2 - 1,
                 'min_rel_margin_A': minA.str(8), 'min_rel_margin_A_at': kA,
                 'min_rel_margin_B': minB.str(8), 'min_rel_margin_B_at': kB,
                 'prec_bits': prec_bits, 'seconds': round(time.time() - t0, 1)})
    return info


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('certificate')
    ap.add_argument('--prec', type=int, default=320)
    ap.add_argument('--out')
    a = ap.parse_args()
    raw = open(a.certificate, 'rb').read()
    d = json.loads(raw)
    info = run(d, a.prec)
    info['certificate_sha256'] = hashlib.sha256(raw).hexdigest()
    s = json.dumps(info, indent=1)
    print(s)
    if a.out:
        open(a.out, 'w').write(s + '\n')
