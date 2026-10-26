"""Exact scores for C1a / C1b step functions (values read as exact rationals), Kronecker substitution with big ints."""
import json, sys
from fractions import Fraction
from math import lcm
mode, path = sys.argv[1], sys.argv[2]
vals = [Fraction(v) for v in json.load(open(path))["values"]]
n = len(vals)
den = lcm(*(v.denominator for v in vals))
a = [int(v * den) for v in vals]
assert min(a) >= 0
def conv(x, y):
    B = (max(x).bit_length() + max(y).bit_length() + len(x).bit_length() + 2)
    X = 0
    for v in reversed(x): X = (X << B) | v
    Y = 0
    for v in reversed(y): Y = (Y << B) | v
    Z = X * Y; mask = (1 << B) - 1
    return [(Z >> (B * k)) & mask for k in range(len(x) + len(y) - 1)]
if mode == "1a":
    c = conv(a, a)
    s = Fraction(2 * n * max(c), sum(a) ** 2)
    print(f"1a n={n} peak_index={c.index(max(c))} score={s}  ~{float(s):.15f}")
else:
    assert sum(vals) == Fraction(n, 2), f"mass off by {sum(vals) - Fraction(n, 2)}"
    assert max(a) <= den
    comp = [den - v for v in reversed(a)]
    c = conv(a, comp)
    s = Fraction(2 * max(c), n * den * den)
    print(f"1b n={n} score~{float(s):.15f}")
