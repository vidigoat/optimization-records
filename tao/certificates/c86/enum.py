"""Enumerate monic integer polynomials of degree d with all roots real in [0, B] and trace <= T.

Robinson-style recursion on derivatives: write P = sum c_i x^i (c_d = 1). For j = 1..d the (d-j)-th
derivative Q_j depends only on c_{d-j}, ..., c_d, and must itself have all roots in [0, B]. Q_j = R_j + c_{d-j} (d-j)!
where R_j is fixed by the earlier coefficients, so the admissible c_{d-j} form an interval, read off from the
values of R_j at 0, B and at the roots of Q_j' = (d-j+1)^{-1}... (the previous level). Intervals are widened by a
small tolerance; every final candidate is re-checked with numpy roots and then exactly with a Sturm count.

    python enum.py d B T  -> prints coefficient lists (low to high) of irreducible candidates
"""
import sys
from fractions import Fraction
from math import factorial

import numpy as np
import sympy

d = int(sys.argv[1])
B = Fraction(sys.argv[2])
T = int(sys.argv[3])
Bf = float(B)
EPS = 1e-7


def deriv_coeffs(c, m):
    """coefficients (low..high) of the m-th derivative of sum c_i x^i, c given low..high."""
    return [c[i] * factorial(i) // factorial(i - m) for i in range(m, len(c))]


def polyval(q, x):
    r = 0.0
    for a in reversed(q):
        r = r * x + a
    return r


out = []


def rec(j, top, crit):
    """top: coefficients c_{d-j+1..d} fixed (list, index 0 = c_{d-j+1}); crit: sorted roots of Q_{j-1} (in [0,B])."""
    if j > d:
        out.append(top)
        return
    m = d - j
    # Q_j = sum_{i>=m} c_i i!/(i-m)! x^{i-m};  constant term c_m * m!
    full = [0] + top                                   # c_m placeholder at index 0
    R = deriv_coeffs([0] * m + full, m)                # low..high of Q_j with c_m = 0
    pts = [0.0] + list(crit) + [Bf]
    sgn_top = 1                                        # leading coefficient positive
    # Q_j has j real roots in [0,B] iff, with critical points r_1<...<r_{j-1} of Q_j (roots of Q_{j-1}),
    # sign Q_j(B) >= 0, sign Q_j(r_{j-1}) <= 0, alternating downwards, and (-1)^j Q_j(0) >= 0.
    lo, hi = -np.inf, np.inf
    for k, x in enumerate(reversed(pts)):              # k = 0 at B
        need = 1 if k % 2 == 0 else -1                 # required sign of Q_j at this point
        v = polyval(R, x)                              # Q_j(x) = v + c_m m!
        f = factorial(m)
        if need == 1:                                  # v + c f >= 0  -> c >= -v/f
            lo = max(lo, -v / f - EPS)
        else:
            hi = min(hi, -v / f + EPS)
    if j == 1:                                         # c_{d-1} = -trace, trace <= T, trace >= 0
        lo, hi = max(lo, -T), min(hi, 0)
    for c in range(int(np.ceil(lo)), int(np.floor(hi)) + 1):
        if j == d and c == 0 and d > 1:
            continue                                   # x | P: reducible unless P = x
        Q = [c * factorial(m)] + R[1:]
        if j == 1:
            roots = [-Q[0] / Q[1]]
        else:
            roots = np.roots(Q[::-1])
            if np.max(np.abs(roots.imag)) > 1e-6:
                continue
            roots = np.sort(roots.real)
            if roots[0] < -1e-7 or roots[-1] > Bf + 1e-7:
                continue
        rec(j + 1, [c] + top, list(np.clip(roots, 0.0, Bf)))


rec(1, [1], [])
x = sympy.symbols("x")
count = 0
for top in out:
    coeffs = top                                       # low..high
    P = sympy.Poly(list(reversed(coeffs)), x)
    if d > 1 and not P.is_irreducible:
        continue
    # exact check: all roots in [0, B] (Sturm counts)
    if P.count_roots(0, B) != d:
        continue
    count += 1
    print(coeffs, flush=True)
print(f"# degree {d}: {count} irreducible polynomials with all roots in [0,{B}] and trace <= {T}", file=sys.stderr)
