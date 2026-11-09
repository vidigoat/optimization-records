"""Exact C1d score of a signed integer step function: 2n * max_m |(v*v)_m| / (sum v)^2.
v = P - Q with P, Q >= 0, so v*v = P*P - 2 P*Q + Q*Q; each product is one big-integer multiplication."""
import json, sys
from fractions import Fraction
d = json.load(open(sys.argv[1])); v = [int(x) for x in d["v_int"]]; n = len(v)
assert n == int(d["n"])
P = [x if x > 0 else 0 for x in v]; Q = [-x if x < 0 else 0 for x in v]
M = max(max(P), max(Q)); W = (2 * M.bit_length() + n.bit_length() + 3 + 7) // 8
pack = lambda a: int.from_bytes(b"".join(x.to_bytes(W, "little") for x in a), "little")
def coeffs(Z):
    b = Z.to_bytes(W * 2 * n, "little"); return [int.from_bytes(b[W * k: W * (k + 1)], "little") for k in range(2 * n - 1)]
XP, XQ = pack(P), pack(Q)
pp, pq, qq = coeffs(XP * XP), coeffs(XP * XQ), coeffs(XQ * XQ)
best, bk = -1, -1
for k in range(2 * n - 1):
    c = abs(pp[k] - 2 * pq[k] + qq[k])
    if c > best: best, bk = c, k
S = sum(v); s = Fraction(2 * n * best, S * S)
print(f"1d n={n} peak_index={bk} |peak|={best} S={S}\nscore={s.numerator}/{s.denominator}\n~{float(s):.15f}")
