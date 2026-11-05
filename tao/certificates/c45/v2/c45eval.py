"""Exact evaluation of the Romanoff obstruction bound from a base histogram (bins nu0 <= C) and cluster primes.

    S = sum_{nu0} (delta0(nu0)/M0) * sum_H P(H) min(1/2, nu0 H kappa),   kappa = M / (T phi(M) lam),
with M = M0 * prod(cluster primes), T = T0 * prod(cluster orders), lam < log 2 rational, and H = prod_d h_d the
cluster factor (h_d = number of r mod d not killed by the order-d cluster primes; exact occupancy law).
Bins above C are given only by their total count; they are valid only if (C+1) * H_min * kappa >= 1/2, which is
checked, and then they contribute exactly 1/2 each.

    python c45eval.py spec.json hist.txt [C]
spec: {"base": [...], "clusters": {"23": [47, 178481], ...}}
"""
import json, sys
from fractions import Fraction
from math import lcm, prod


def ord2(q):
    x, d = 2 % q, 1
    while x != 1:
        x, d = x * 2 % q, d + 1
    return d


def lam_log2():   # 2 sum_{j<60} 1/((2j+1) 3^(2j+1)) < log 2
    return 2 * sum(Fraction(1, (2 * j + 1) * 3 ** (2 * j + 1)) for j in range(60))


def cluster_law(d, qs):
    """Exact law of h = d - #distinct killed classes, as {h: count}, counts summing to prod(qs)."""
    st = {0: 1}                                  # number of distinct killed classes -> count
    for q in qs:
        nx = {}
        for j, c in st.items():
            nx[j] = nx.get(j, 0) + c * (q - d + j)          # kills nothing, or a class already killed
            if j < d: nx[j + 1] = nx.get(j + 1, 0) + c * (d - j)
        st = nx
    return {d - j: c for j, c in st.items()}


def main():
    spec = json.load(open(sys.argv[1])); hist = sys.argv[2]
    base = spec["base"]; cl = {int(k): v for k, v in spec["clusters"].items()}
    for d, qs in cl.items():
        assert all(ord2(q) == d for q in qs), d
    M0 = prod(base); T0 = lcm(*[ord2(q) for q in base])
    for d in cl: assert T0 % d != 0
    Mc = prod(q for qs in cl.values() for q in qs); M = M0 * Mc
    T = T0 * prod(cl); phi = prod(q - 1 for q in base) * prod(q - 1 for qs in cl.values() for q in qs)
    lam = lam_log2(); kappa = Fraction(M, T * phi) / lam
    # combined cluster law
    Hl = {1: 1}
    for d, qs in sorted(cl.items()):
        law = cluster_law(d, qs); nx = {}
        for h1, c1 in Hl.items():
            for h2, c2 in law.items():
                nx[h1 * h2] = nx.get(h1 * h2, 0) + c1 * c2
        Hl = nx
    assert sum(Hl.values()) == Mc
    hs = sorted(Hl); Hmin = hs[0]
    # bins
    bins = []; above = None
    for line in open(hist):
        if line.startswith('#'):
            continue
        n, c = line.split(); bins.append((int(n), int(c)))
    tot = sum(c for _, c in bins)
    above = M0 - tot
    assert above >= 0
    Cmax = max(n for n, _ in bins)
    if above > 0:
        C = int(sys.argv[3]) if len(sys.argv) > 3 else Cmax
        assert Cmax <= C
        assert (C + 1) * Hmin * kappa >= Fraction(1, 2), "cut C too small for the regime 1/2"
    # S = sum over bins of delta0/M0 * g(nu0) + above/(2 M0); g grouped by J(nu0) = #{h_j : 2 nu0 kappa h_j < 1}
    pre_ph = [0]; pre_p = [0]
    for h in hs:
        pre_ph.append(pre_ph[-1] + Hl[h] * h); pre_p.append(pre_p[-1] + Hl[h])
    Kn, Kd = kappa.numerator, kappa.denominator
    J = len(hs); acc = {}                        # J -> [sum delta0*nu0, sum delta0]
    for n, c in sorted(bins):
        while J > 0 and not (2 * n * Kn * hs[J - 1] < Kd):
            J -= 1
        a = acc.setdefault(J, [0, 0]); a[0] += c * n; a[1] += c
    S = Fraction(above, 2 * M0)
    for Jv, (sn, sc) in acc.items():
        S += (kappa * sn * Fraction(pre_ph[Jv], Mc) + Fraction(sc * (Mc - pre_p[Jv]), 2 * Mc)) / M0
    print(f"base primes {len(base)}, cluster primes {sum(len(v) for v in cl.values())}, T0 {T0}, bins {len(bins)}, above-C count {'0' if above == 0 else 'yes'}")
    print(f"C_min for the cut: {(Fraction(1, 2) / (Hmin * kappa)).__ceil__() - 1}")
    from decimal import Decimal, getcontext
    getcontext().prec = 60
    print("S =", Decimal(S.numerator) / Decimal(S.denominator))
    json.dump({"S_num": str(S.numerator), "S_den": str(S.denominator)}, open(hist + ".S.json", "w"))


if __name__ == "__main__":
    main()
