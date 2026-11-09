"""Exact C1a score of an integer step function: 2n * max_k (a*a)_k / (sum a)^2, via one big-integer product."""
import json, sys
from fractions import Fraction
d = json.load(open(sys.argv[1])); a = [int(x) for x in d["v_int"]]; n = len(a)
assert n == int(d["n"]) and min(a) >= 0
W = (2 * max(a).bit_length() + n.bit_length() + 2 + 7) // 8  # bytes per slot
X = int.from_bytes(b"".join(v.to_bytes(W, "little") for v in a), "little")
Zb = (X * X).to_bytes(W * (2 * n), "little")
best, bk = -1, -1
for k in range(2 * n - 1):
    c = int.from_bytes(Zb[W * k: W * (k + 1)], "little")
    if c > best: best, bk = c, k
S = sum(a); s = Fraction(2 * n * best, S * S)
print(f"1a n={n} peak_index={bk} peak={best} S={S}\nscore={s.numerator}/{s.denominator}\n~{float(s):.15f}")
