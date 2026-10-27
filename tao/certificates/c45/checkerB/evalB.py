#!/usr/bin/env python3
"""evalB.py -- independent exact evaluator ("checker B") for the bound

    dbar <= S = sum_nu delta(nu) * min( 1/(2M), nu / (T * phi(M) * lambda) ),

where P = B u C_1 u ... u C_k, M = prod P, T = lcm of the orders ord_q(2),
nu(a) = #{ r mod T : a mod q != 2^r mod q for all q in P }, delta(nu) = #{a mod M :
nu(a) = nu}, and lambda is a rational number with lambda < log 2 (replacing log 2
by a smaller number can only increase S, so S is still an upper bound).

Structure used (all of it is verified below):
  * every base prime has order dividing T0 = lcm of the base orders;
  * every cluster C_j consists of primes with the same prime order p_j,
    gcd(p_j, T0) = 1 and the p_j are pairwise distinct.
Then T = T0 * prod p_j, and by CRT nu(a) = nu0(a mod M0) * prod_j h_j(a mod M_Cj),
where h_j(b) = #{ r mod p_j : b mod q != 2^r mod q for all q in C_j }.  The residues
a mod M0, a mod M_C1, ... are independent, so delta is the multiplicative
convolution of the base histogram delta0 (read from a file produced by histB)
and the cluster distributions delta_j.

Cluster distribution (occupancy formula).  For q in C (order p), b mod q is either
a non-power of 2 (q - p residues, removes nothing) or 2^s for exactly one class
s mod p.  If the set A of primes whose residue is a power of 2 has |A| = m, the
number of ways the m chosen classes cover exactly k distinct classes is
binom(p, k) * k! * S2(m, k) (S2 = Stirling numbers of the second kind).  Hence
    delta_C(p - k) = sum_m e_m * binom(p, k) * k! * S2(m, k),
    e_m = sum_{|A| = m} prod_{q not in A} (q - p)   (coefficients of prod (x + q - p)).
This is cross-checked by (i) an independent DP over the number of occupied
classes and (ii) literal enumeration of b mod M_C for small clusters.

lambda: partial sum of log 2 = -log(1 - 1/2) = sum_{k>=1} 1/(k 2^k).  All terms are
positive, so every partial sum is < log 2, and the omitted tail is
sum_{k>N} 1/(k 2^k) < (1/(N+1)) * sum_{k>N} 2^-k = 1/((N+1) 2^N).  With N = 200
the error is < 4e-63.  (This deliberately differs from the atanh series.)

Usage: evalB.py spec.json [--brute-limit=N]
   spec.json = {"base": [...], "hist": "path", "clusters": [[...], [...], ...]}
   a relative "hist" path is resolved against the spec file's directory first,
   then against the current directory.
Python 3 standard library only.
"""
import sys, json, os, math, bisect, time
from fractions import Fraction

# ------------------------------------------------------------------ arithmetic helpers

def is_prime(n):
    if n < 2:
        return False
    small = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
    for p in small:
        if n % p == 0:
            return n == p
    if n < 10**12:
        i = 41
        while i * i <= n:
            if n % i == 0:
                return False
            i += 2
        return True
    # deterministic Miller-Rabin for n < 3.3e24
    if n >= 3317044064679887385961981:
        raise ValueError("number too large for deterministic primality test")
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2; s += 1
    for a in small:
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True

def factorize(n):
    f = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f

def order_of_2(q):
    """ord_q(2): smallest k > 0 with 2^k = 1 mod q, from the factorisation of q-1,
    certified by 2^k = 1 and 2^(k/r) != 1 for every prime r | k.  For q < 2e7 it
    is additionally confirmed by direct iteration of x -> 2x mod q."""
    k = q - 1
    for r in factorize(q - 1):
        while k % r == 0 and pow(2, k // r, q) == 1:
            k //= r
    assert pow(2, k, q) == 1
    for r in factorize(k):
        assert pow(2, k // r, q) != 1
    if q < 2 * 10**7:
        x, j = 2 % q, 1
        while x != 1:
            x = x * 2 % q
            j += 1
        assert j == k, (q, j, k)
    return k

def lcm(a, b):
    return a // math.gcd(a, b) * b

def log2_lower(N=200):
    """Rational lambda = sum_{k=1}^N 1/(k 2^k) < log 2, error < 1/((N+1) 2^N)."""
    s = Fraction(0)
    for k in range(1, N + 1):
        s += Fraction(1, k * 2**k)
    return s, Fraction(1, (N + 1) * 2**N)

# ------------------------------------------------------------------ cluster distributions

def stirling2_table(mmax):
    S = [[0] * (mmax + 1) for _ in range(mmax + 1)]
    S[0][0] = 1
    for m in range(1, mmax + 1):
        for k in range(1, m + 1):
            S[m][k] = k * S[m - 1][k] + S[m - 1][k - 1]
    return S

def cluster_dist_occupancy(primes, p):
    """{h: count} via the occupancy / Stirling formula."""
    m = len(primes)
    # e[j] = coefficient of x^j in prod (x + (q - p))
    e = [1]
    for q in primes:
        w = q - p
        ne = [0] * (len(e) + 1)
        for j, c in enumerate(e):
            ne[j] += c * w
            ne[j + 1] += c
        e = ne
    S = stirling2_table(m)
    dist = {}
    for j in range(m + 1):
        for k in range(0, min(j, p) + 1):
            if S[j][k] == 0:
                continue
            cnt = e[j] * math.comb(p, k) * math.factorial(k) * S[j][k]
            dist[p - k] = dist.get(p - k, 0) + cnt
    return dist

def cluster_dist_dp(primes, p):
    """Independent check: DP over primes with state = number of occupied classes."""
    st = {0: 1}
    for q in primes:
        ns = {}
        for k, c in st.items():
            ns[k] = ns.get(k, 0) + c * (q - p)       # not a power of 2
            if k:
                ns[k] = ns.get(k, 0) + c * k          # hits an occupied class
            if k < p:
                ns[k + 1] = ns.get(k + 1, 0) + c * (p - k)  # hits a new class
        st = ns
    return {p - k: c for k, c in st.items()}

def cluster_dist_brute(primes, p):
    """Literal enumeration of b mod M_C (M_C small): h(b) = #{r mod p : b mod q != 2^r mod q}."""
    MC = 1
    for q in primes:
        MC *= q
    # residue maps computed literally from the definition (r over Z_p, 2^r mod q)
    pw = [{pow(2, r, q): r for r in range(p)} for q in primes]
    for q, d in zip(primes, pw):
        assert len(d) == p  # 2^r mod q distinct for r mod p (p = order)
    dist = {}
    for b in range(MC):
        hit = set()
        for q, d in zip(primes, pw):
            r = d.get(b % q)
            if r is not None:
                hit.add(r)
        h = p - len(hit)
        dist[h] = dist.get(h, 0) + 1
    return dist

def cluster_dist_grouped(primes, p):
    """Semi-literal check: for every prime, group the residues b mod q by the class
    {r mod p : 2^r = b mod q} they remove (computed literally), then enumerate all
    combinations of groups and take the size of the union directly."""
    import itertools
    groups = []
    for q in primes:
        g = {}
        hitr = {}
        for r in range(p):
            hitr[pow(2, r, q)] = r
        for b in range(q) if q <= 10**6 else []:
            key = hitr.get(b)
            g[key] = g.get(key, 0) + 1
        if q > 10**6:  # same grouping without scanning all residues
            g = {r: 1 for r in hitr.values()}
            g[None] = q - len(hitr)
        groups.append(list(g.items()))
    dist = {}
    for combo in itertools.product(*groups):
        cnt = 1
        hit = set()
        for key, c in combo:
            cnt *= c
            if key is not None:
                hit.add(key)
        h = p - len(hit)
        dist[h] = dist.get(h, 0) + cnt
    return dist

# ------------------------------------------------------------------ main evaluation

def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    opts = dict(a[2:].split("=", 1) if "=" in a else (a[2:], "1") for a in sys.argv[1:] if a.startswith("--"))
    brute_limit = int(opts.get("brute-limit", 3 * 10**6))
    if len(args) != 1:
        print(__doc__)
        sys.exit(2)
    t_start = time.time()
    spec_path = args[0]
    spec = json.load(open(spec_path))
    base = [int(x) for x in spec["base"]]
    clusters = [[int(x) for x in c] for c in spec["clusters"]]
    hist_path = spec["hist"]
    if not os.path.isabs(hist_path):
        cand = os.path.join(os.path.dirname(os.path.abspath(spec_path)), hist_path)
        if os.path.exists(cand):
            hist_path = cand
    print("histogram file   :", os.path.abspath(hist_path))

    # ---- primes, distinctness, orders
    allp = base + [q for c in clusters for q in c]
    assert len(set(allp)) == len(allp), "primes are not distinct"
    for q in allp:
        assert q % 2 == 1 and is_prime(q), "not an odd prime: %d" % q
    ords = {q: order_of_2(q) for q in allp}
    T0 = 1
    for q in base:
        T0 = lcm(T0, ords[q])
    M0 = 1
    phiM0 = 1
    for q in base:
        M0 *= q
        phiM0 *= q - 1
    print("base primes      :", base)
    print("base orders      :", [ords[q] for q in base])
    print("T0 = %d, M0 = %d" % (T0, M0))
    cl_p = []
    for c in clusters:
        os_ = {ords[q] for q in c}
        assert len(os_) == 1, "cluster %s has different orders %s" % (c, os_)
        p = os_.pop()
        assert is_prime(p), "cluster order %d is not prime" % p
        assert math.gcd(p, T0) == 1, "cluster order %d not coprime to T0" % p
        cl_p.append(p)
    assert len(set(cl_p)) == len(cl_p), "cluster orders not distinct"
    print("cluster orders   :", cl_p)
    T = T0
    for p in cl_p:
        T *= p
    Tl = 1
    for q in allp:
        Tl = lcm(Tl, ords[q])
    assert T == Tl, "T != lcm of all orders"
    M = 1
    phiM = 1
    for q in allp:
        M *= q
        phiM *= q - 1
    print("T = %d" % T)
    print("M = %d (%d digits), phi(M) has %d digits" % (M, len(str(M)), len(str(phiM))))

    # ---- base histogram
    hist = []
    prev = -1
    with open(hist_path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            a, b = line.split()
            nu, dl = int(a), int(b)
            assert nu > prev, "histogram not strictly ascending"
            assert 0 <= nu <= T0 and dl > 0
            prev = nu
            hist.append((nu, dl))
    s0 = sum(d for _, d in hist)
    s1 = sum(d * n for n, d in hist)
    s2 = sum(d * n * n for n, d in hist)
    assert s0 == M0, "sum delta0 != M0"
    assert s1 == T0 * phiM0, "sum delta0*nu != T0*phi(M0)"
    # second moment: sum_a nu^2 = T0 * sum_{x in Z_T0} prod_q (q - 2 + [d_q | x])
    e2 = 0
    for x in range(T0):
        t = 1
        for q in base:
            t *= q - 2 + (1 if x % ords[q] == 0 else 0)
        e2 += t
    e2 *= T0
    assert s2 == e2, "second moment check failed"
    print("base histogram   : %d entries, sum = M0 OK, first moment OK, second moment OK" % len(hist))

    # ---- cluster distributions
    cdists = []
    for c, p in zip(clusters, cl_p):
        d1 = cluster_dist_occupancy(c, p)
        d2 = cluster_dist_dp(c, p)
        assert d1 == d2, "occupancy formula and DP disagree for %s" % c
        MC = 1
        phiC = 1
        for q in c:
            MC *= q
            phiC *= q - 1
        assert sum(d1.values()) == MC
        assert sum(h * v for h, v in d1.items()) == p * phiC
        note = "formula = DP"
        ncomb = 1
        for q in c:
            ncomb *= p + 1
        if ncomb <= 2 * 10**6:
            d4 = cluster_dist_grouped(c, p)
            assert d1 == d4, "grouped enumeration disagrees for %s" % c
            note += " = grouped enumeration (%d combos)" % ncomb
        if MC <= brute_limit:
            d3 = cluster_dist_brute(c, p)
            assert d1 == d3, "brute force disagrees for %s" % c
            note += " = brute force over %d residues" % MC
        print("cluster %-28s p=%-4d h-dist %s  [%s]" % (c, p, dict(sorted(d1.items(), reverse=True)), note))
        cdists.append(d1)

    # ---- multiplier distribution D(H), H = prod h_j
    D = {1: 1}
    for d in cdists:
        nd = {}
        for H, c in D.items():
            for h, v in d.items():
                k = H * h
                nd[k] = nd.get(k, 0) + c * v
        D = nd
    MCall = 1
    for c in clusters:
        for q in c:
            MCall *= q
    assert sum(D.values()) == MCall
    Hs = sorted(D)
    print("multiplier distribution: %d distinct products" % len(Hs))
    PD = [0] * (len(Hs) + 1)
    PH = [0] * (len(Hs) + 1)
    for i, H in enumerate(Hs):
        PD[i + 1] = PD[i] + D[H]
        PH[i + 1] = PH[i] + D[H] * H
    TotD = PD[-1]
    assert TotD * M0 == M
    # global first-moment check: sum_a nu(a) = T * phi(M)
    assert s1 * PH[-1] == T * phiM

    def evaluate(lam):
        # term(nu, H) = min(1/(2M), nu H / (T phi lam)); nu H/(T phi lam) < 1/(2M)
        # <=> 2 M nu H * lam.den < T phi lam.num
        X = T * phiM * lam.numerator
        Y = 2 * M * lam.denominator
        A = 0  # sum delta0 * nu * sum_{H below threshold} D(H) H
        Bc = 0  # sum delta0 * sum_{H at/above threshold} D(H)
        for nu, dl in hist:
            if nu == 0:
                continue
            den = Y * nu
            t = -(-X // den)  # H < X/den  <=>  H < ceil(X/den) for integer H
            k = bisect.bisect_left(Hs, t)
            A += dl * nu * PH[k]
            Bc += dl * (TotD - PD[k])
        return Fraction(A) / (T * phiM * lam) + Fraction(Bc, 2 * M)

    lam, err = log2_lower(200)
    assert err < Fraction(1, 10**40)
    S = evaluate(lam)
    S_hi_lam = evaluate(lam + err)  # value with a lambda >= log 2 (for information only)
    num, den = S.numerator, S.denominator
    print()
    print("lambda = sum_{k=1}^{200} 1/(k 2^k)  (< log 2, error < %.1e)" % float(err))
    print("S (exact fraction) = %d / %d" % (num, den))
    # 40 decimal digits (truncated)
    digits = 40
    q = num * 10**digits // den
    s = str(q).rjust(digits + 1, "0")
    print("S = %s.%s...  (truncated to %d decimals)" % (s[:-digits], s[-digits:], digits))
    up = -(-num * 10**15 // den)
    su = str(up).rjust(16, "0")
    print("S rounded UP to 15 decimals: %s.%s" % (su[:-15], su[-15:]))
    diff = S - S_hi_lam
    print("sensitivity: S(lambda) - S(lambda + err) = %.3e (>= 0 expected)" % float(diff))
    assert diff >= 0
    print("elapsed %.1fs" % (time.time() - t_start))

if __name__ == "__main__":
    main()
