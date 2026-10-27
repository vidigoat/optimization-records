"""Checker A, exact evaluation stage.

Given a base prime set with its exact histogram delta0(nu0) (nu0 counted modulo T0 = lcm of base orders)
and clusters of primes whose orders are pairwise coprime and coprime to T0, evaluate exactly

  S = sum_{a mod M} min( 1/(2M), nu(a) / (T phi(M) lambda) ),   lambda <= log 2 rational,

where by CRT nu(a) = nu0(a mod M0) * prod_c h_c(a mod M_c), T = T0 * prod_c T_c.
All arithmetic is exact integer / Fraction arithmetic.
"""
from fractions import Fraction
from math import gcd, prod
from bisect import bisect_right
from collections import Counter
import sys, json


def ord2(q):
    x, d = 2 % q, 1
    while x != 1:
        x = x * 2 % q
        d += 1
    return d


def lcm(a, b):
    return a // gcd(a, b) * b


def log2_lower(terms=60):
    # log 2 = 2 * sum_{j>=0} 1/((2j+1) 3^(2j+1)); positive terms, so partial sums are lower bounds
    return sum(Fraction(2, (2 * j + 1) * 3 ** (2 * j + 1)) for j in range(terms))


def cluster_dist(primes):
    """Distribution of h = #{r mod d : r survives all primes of the cluster}; all primes have order d.
    For each prime q: with weight q-d the residue is not a power of 2 (no removal); otherwise it removes
    one class r = s mod d, each s with weight 1. Track the number k of distinct removed classes."""
    d = ord2(primes[0])
    dist = {0: 1}
    for q in primes:
        assert ord2(q) == d
        nd = Counter()
        for k, c in dist.items():
            nd[k] += c * (q - d)       # not a power of two
            nd[k] += c * k             # removes an already removed class
            if k < d:
                nd[k + 1] += c * (d - k)  # removes a new class
        dist = dict(nd)
    out = {d - k: c for k, c in dist.items() if c}
    assert sum(out.values()) == prod(primes)
    return d, out


def load_hist(path):
    h = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            a, b = line.split()
            a, b = int(a), int(b)
            assert a not in h and b > 0
            h[a] = b
    return h


def evaluate(base_primes, hist, clusters, log_terms=60, verbose=True):
    base_primes = list(base_primes)
    assert len(set(base_primes)) == len(base_primes)
    T0 = 1
    for q in base_primes:
        T0 = lcm(T0, ord2(q))
    M0 = prod(base_primes)
    assert sum(hist.values()) == M0, "histogram total != M0"
    assert all(0 <= nu <= T0 for nu in hist)
    # clusters
    mult = Counter({1: 1})
    T = T0
    allp = list(base_primes)
    for c in clusters:
        d, dist = cluster_dist(c)
        assert gcd(d, T) == 1, "cluster order not coprime to the rest"
        T *= d
        allp += list(c)
        new = Counter()
        for m, w in mult.items():
            for h, c2 in dist.items():
                new[m * h] += w * c2
        mult = new
    assert len(set(allp)) == len(allp)
    M = prod(allp)
    phi = prod(q - 1 for q in allp)
    T_check = 1
    for q in allp:
        T_check = lcm(T_check, ord2(q))
    assert T_check == T
    assert sum(mult.values()) * M0 == M
    lam = log2_lower(log_terms)
    A = Fraction(1, 2 * M)
    B = Fraction(1, 1) / (T * phi * lam)    # value per unit nu: nu * B
    ms = sorted(mult)
    pw = [0]
    pmw = [0]
    for m in ms:
        pw.append(pw[-1] + mult[m])
        pmw.append(pmw[-1] + mult[m] * m)
    W = pw[-1]
    S = Fraction(0)
    lin = 0     # integer accumulators: S = A*cnt_const + B*lin
    cnt_const = 0
    for nu, dlt in hist.items():
        if nu == 0:
            continue
        # term uses linear branch iff nu*m*B <= A  <=>  m <= A/(B nu)
        thr = A / (B * nu)
        k = bisect_right(ms, thr.numerator // thr.denominator)
        lin += dlt * nu * pmw[k]
        cnt_const += dlt * (W - pw[k])
    S = A * cnt_const + B * lin
    info = dict(n_primes=len(allp), M=M, phi=phi, T=T, T0=T0, base_entries=len(hist), mult_entries=len(ms))
    return S, info


def ceil_dec(fr, places):
    s = 10 ** places
    return -((-fr.numerator * s) // fr.denominator)


if __name__ == '__main__':
    spec = json.load(open(sys.argv[1]))
    hist = load_hist(spec['hist'])
    S, info = evaluate(spec['base'], hist, spec['clusters'])
    print(json.dumps({k: str(v) for k, v in info.items()}, indent=1))
    print('S =', S.numerator, '/', S.denominator)
    up = ceil_dec(S, 15)
    print('S (decimal, 30 digits) ~ ', S.numerator * 10 ** 30 // S.denominator)
    print('proved upper bound (rounded up, 15 places): 0.%015d' % up if up < 10 ** 15 else up)
