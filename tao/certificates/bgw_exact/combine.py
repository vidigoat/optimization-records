"""Combine exact.c outputs into a rigorous lower bound for C16 (mpmath interval arithmetic at the end)."""
import sys
from fractions import Fraction
from mpmath import iv, mp
cnt = s0 = s1 = s2 = 0
for line in open(sys.argv[1]):
    c, a, b, d = map(int, line.split()); cnt += c; s0 += a; s1 += b; s2 += d
ONE = 2 ** 100
# total over all lattice points = 4 x quadrant; each term's true value is in [floor, floor + 1] / 2^100
lo = lambda s: Fraction(4 * s, ONE)
hi = lambda s: Fraction(4 * (s + cnt), ONE)
iv.dps = 50
I = lambda f_lo, f_hi: iv.mpf([iv.mpf(f_lo.numerator) / f_lo.denominator, iv.mpf(f_hi.numerator) / f_hi.denominator])
S0, S1, S2 = I(lo(s0), hi(s0)), I(lo(s1), hi(s1)), I(lo(s2), hi(s2))
delta = S2 / S1
F = S0 * S0 / (iv.pi * S1) - iv.log(delta) - iv.log(1 + iv.log(delta))
mp.dps = 25
print(f"lattice points: {4 * cnt}")
lo_, hi_ = lambda x: mp.mpf(x.a.a), lambda x: mp.mpf(x.b.b)
print(f"delta in [{mp.nstr(lo_(delta), 15)}, {mp.nstr(hi_(delta), 15)}]  (>= 1 required)")
print(f"C16 >= {mp.nstr(lo_(F), 18)}   (interval width {mp.nstr(hi_(F) - lo_(F), 3)})")
