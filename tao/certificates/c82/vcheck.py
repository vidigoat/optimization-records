"""Rigorous Arb upper bound for the Doche limit-point value B of a C82 certificate (Pinho's JSON format).

Written for this project from the mathematics (it imports only json/fractions/multiprocessing/flint).

    B = (1/D) [ sum_j m_j log M(R_j o T) + int_0^1 max(0, log|A(T(u))| - log|Q(T(u))|) dt ],   u = exp(2 pi i t),

T(z) = z(1-z), Q = prod R_j^{m_j}, A = prod S_i^{q_i}, D = deg Q = sum m_j deg R_j.

Exact checks: integer coefficients, |lead R_j| = 1, m_j positive integers, q_i >= 0 rationals, deg A < D,
gcd(S_i, R_j) = 1 whenever q_i > 0 (FLINT exact gcd).

Every polynomial P = R_j o T or S_i o T is composed exactly in Z[z]. Arb (FLINT complex_roots) encloses its complete
root multiset; the multiplicities are checked to add up to deg P. Hence for every complex u
    log|P(u)| = log|lead P| + sum_r log|u - r|      (r over the enclosed roots, with multiplicity).

Mahler part (Jensen): log M(P) = log|lead P| + sum_r log max(1, |r|).

Integral part (only an upper bound is needed).  g(t) = sum_r c_r log|u(t) - r| + const, with c_r = q_i/D for roots
of S_i o T and c_r = -m_j/D for roots of R_j o T.  g(1-t) = g(t) because all coefficients are real, so the
integral over [0,1] is twice the integral over [0,1/2].  [0,1/2] is tiled exactly by dyadic panels with rational
end points.  For a panel of width w and midpoint c, |u(t) - u(c)| <= 2 pi |t - c| <= pi w, so with
d_r = |u(c) - r| every distance on the panel lies in [d_r - pi w, d_r + pi w].

  (1) g_up = const + sum_{c_r>0} c_r log(d_r + pi w) + sum_{c_r<0} c_r log(d_r - pi w)   (needs d_r - pi w > 0 for
      c_r < 0).  If g_up < 0 the panel contributes 0.
  (2) Otherwise, if d_r - pi w > 0 for all roots, g is smooth on the panel with
          g'(t) = Re sum c_r 2 pi i u/(u - r),     g''(t) = Re sum c_r (2 pi i)^2 (-r u)/(u - r)^2,
      so |g''| <= Mb := 4 pi^2 sum |c_r| |r| / (d_r - pi w)^2.  With h = w/2 and reach = |g'(c)| h + Mb h^2/2:
        g(c) - reach > 0:  integral <= w g(c) + Mb w^3/24        (Taylor, linear term integrates to 0)
        g(c) + reach < 0:  integral = 0
        otherwise:         integral <= w max(0,g(c)) + |g'(c)| w^2/4 + Mb w^3/24   (positive part is 1-Lipschitz)
      The panel is accepted when the Taylor slack is below tol*w; otherwise (or if (2) does not apply) it is bisected.
All quantities are Arb balls; panel results are upper end points passed back exactly (mantissa, exponent).

The 2^level top-level panels are processed in parallel (the exact dyadic tiling is unchanged).

usage: python vcheck.py certificate.json [--tol 1e-10] [--prec 128] [--target X] [--procs N] [--out result.json]
"""
import argparse, json, time, hashlib
from fractions import Fraction
from multiprocessing import Pool
from flint import arb, acb, fmpz_poly, ctx


def frac_arb(q):
    q = Fraction(q)
    return arb(q.numerator) / arb(q.denominator)


def compose_T(coeffs_desc):
    T = fmpz_poly([0, 1, -1])
    p = fmpz_poly([0])
    for c in coeffs_desc:
        p = p * T + int(c)
    return p


def all_roots(P):
    rts = P.complex_roots()
    assert sum(m for _, m in rts) == P.degree(), 'root enclosure incomplete'
    return [r for r, m in rts for _ in range(m)]


def load_and_check(path):
    raw = open(path, 'rb').read()
    cert = json.loads(raw)
    den, num = [], []
    for r in cert['denominator']:
        cs, m = r['coeffs'], r['multiplicity']
        assert all(type(x) is int for x in cs) and len(cs) >= 2, 'bad denominator coeffs'
        assert abs(cs[0]) == 1, 'denominator not monic up to sign'
        assert type(m) is int and m > 0, 'bad multiplicity'
        den.append((r.get('name', '?'), cs, m))
    for s in cert['numerator']:
        cs, q = s['coeffs'], Fraction(s['q_rational'])
        assert all(type(x) is int for x in cs) and len(cs) >= 2 and cs[0] != 0, 'bad numerator coeffs'
        assert q >= 0, 'negative exponent'
        if q > 0:
            num.append((s.get('name', '?'), cs, q))
    D = sum((len(cs) - 1) * m for _, cs, m in den)
    degA = sum((len(cs) - 1) * q for _, cs, q in num)
    assert D > 0 and degA < D, 'degree condition deg A < D fails'
    for _, cs, q in num:
        P = fmpz_poly(cs[::-1])
        for _, cr, _ in den:
            assert P.gcd(fmpz_poly(cr[::-1])).degree() == 0, 'numerator and denominator not coprime'
    return cert, den, num, D, degA, hashlib.sha256(raw).hexdigest()


G = {}


def setup(path, prec, tol):
    ctx.prec = prec
    cert, den, num, D, degA, sha = load_and_check(path)
    Darb = arb(D)
    groups = []      # (signed weight, [roots])
    const = arb(0)
    for _, cs, q in num:
        P = compose_T(cs); wgt = frac_arb(q) / Darb
        const += wgt * arb(abs(int(P.coeffs()[-1]))).log()
        groups.append((wgt, all_roots(P)))
    for _, cs, m in den:
        P = compose_T(cs); wgt = -arb(m) / Darb
        const += wgt * arb(abs(int(P.coeffs()[-1]))).log()
        groups.append((wgt, all_roots(P)))
    G.update(groups=groups, const=const, tol=frac_arb(Fraction(tol)), pi=arb.pi(), maxdepth=60,
             absr=[[abs(r) for r in rs] for _, rs in groups])
    return cert, den, num, D, degA, sha


def mahler(den, D):
    tot = arb(0)
    for _, cs, m in den:
        P = compose_T(cs)
        v = arb(abs(int(P.coeffs()[-1]))).log()
        for r in all_roots(P):
            v += abs(r).max(arb(1)).log()
        tot += arb(m) * v
    return tot / arb(D)


def panel_integral(k, depth, stats):
    """Arb upper bound for the integral over [k/2^depth, (k+1)/2^depth] of max(0, g)."""
    pi, groups, const = G['pi'], G['groups'], G['const']
    w = arb(1) / arb(2) ** depth
    c = arb(2 * k + 1) / arb(2) ** (depth + 1)
    ang = 2 * pi * c
    u = acb(ang.cos(), ang.sin())
    rad = pi * w
    # (1) whole-panel upper bound of g
    gup = const
    separated = True
    dists = []
    for wgt, rts in groups:
        ds = [abs(u - r) for r in rts]
        dists.append(ds)
        if wgt > 0:
            prod = arb(1)
            for d in ds:
                prod *= d + rad
                if not (d - rad > 0):
                    separated = False
            gup += wgt * prod.log()
        else:
            prod = arb(1)
            for d in ds:
                lo = d - rad
                if not (lo > 0):
                    gup = None
                    break
                prod *= lo
            if gup is None:
                separated = False
                break
            gup += wgt * prod.log()
    if gup is not None and gup < 0:
        stats[0] += 1
        stats[4] = max(stats[4], depth)
        return arb(0)
    if separated:
        gc = const; gp = acb(0); Mb = arb(0)
        for (wgt, rts), ds, ar in zip(groups, dists, G['absr']):
            prod = arb(1); s1 = acb(0); s2 = arb(0)
            for r, d, a in zip(rts, ds, ar):
                prod *= d
                s1 += 1 / (u - r)
                lo = d - rad
                s2 += a / (lo * lo)
            gc += wgt * prod.log()
            gp += wgt * s1
            Mb += abs(wgt) * s2
        gpa = abs((acb(0, 2 * pi) * u * gp).real)       # |g'(c)|
        Mb = Mb * 4 * pi * pi                            # >= sup |g''| on the panel
        h = w / 2
        reach = gpa * h + Mb * h * h / 2
        if gc - reach > 0:
            err = Mb * w ** 3 / 24; val = gc * w + err; kind = 1
        elif gc + reach < 0:
            err = arb(0); val = arb(0); kind = 0
        else:
            err = gpa * w * w / 4 + Mb * w ** 3 / 24; val = gc.max(arb(0)) * w + err; kind = 2
        if err < G['tol'] * w:
            stats[kind] += 1
            stats[4] = max(stats[4], depth)
            return val
    assert depth < G['maxdepth'], f'unresolved panel {k}/2^{depth}'
    stats[3] += 1
    return panel_integral(2 * k, depth + 1, stats) + panel_integral(2 * k + 1, depth + 1, stats)


def job(args):
    k, depth = args
    stats = [0, 0, 0, 0, 0]
    t0 = time.time()
    v = panel_integral(k, depth, stats)
    m, e = v.upper().man_exp()        # exact upper end point
    return k, int(m), int(e), stats, time.time() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cert')
    ap.add_argument('--tol', default='1e-10')
    ap.add_argument('--prec', type=int, default=128)
    ap.add_argument('--target', default=None)
    ap.add_argument('--procs', type=int, default=8)
    ap.add_argument('--level', type=int, default=8, help='2^level top-level panels tile [0,1/2]')
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    t0 = time.time()
    cert, den, num, D, degA, sha = setup(a.cert, a.prec, a.tol)
    M = mahler(den, D)
    print('D', D, 'deg A / D', float(degA / D), 'numerators', len(num), 'denominators', len(den), flush=True)
    print('normalised Mahler part', M, flush=True)
    depth = a.level + 1
    tasks = [(k, depth) for k in range(2 ** a.level)]
    with Pool(a.procs, initializer=setup, initargs=(a.cert, a.prec, a.tol)) as pool:
        res = pool.map(job, tasks, chunksize=1)
    ctx.prec = a.prec
    assert sorted(k for k, *_ in res) == list(range(2 ** a.level))
    total = arb(0)
    stats = [0, 0, 0, 0, 0]
    slow = max(r[4] for r in res)
    for k, m, e, st, sec in res:
        total += arb(m) * arb(2) ** e        # exact
        for i in range(4):
            stats[i] += st[i]
        stats[4] = max(stats[4], st[4])
    integral_up = (2 * total).upper()
    B = M + integral_up
    B_up = B.upper()
    res = dict(certificate=a.cert, sha256=sha, D=D, degA_over_D=float(degA / D), numerators_active=len(num),
               mahler_normalized=M.str(25), integral_upper=integral_up.str(25, radius=False),
               B_ball=B.str(25), B_upper=B_up.str(25, radius=False), tol=a.tol, prec=a.prec,
               panels=dict(zero=stats[0], positive=stats[1], crossing=stats[2], bisections=stats[3], max_depth=stats[4]),
               slowest_top_panel_s=round(slow, 1), seconds=round(time.time() - t0, 1))
    if a.target is not None:
        res['target'] = a.target
        res['below_target'] = bool(B_up < frac_arb(Fraction(a.target)))
    print(json.dumps(res, indent=1), flush=True)
    if a.out:
        json.dump(res, open(a.out, 'w'), indent=1)
    if a.target is not None and not res['below_target']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
