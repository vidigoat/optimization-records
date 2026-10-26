"""Second, independent exact C1c scorer: Kronecker substitution with Python big integers.
X = sum w_i 2^(B i), Y = sum w_i 2^(B (n-1-i)); the coefficient of 2^(B (n-1+k)) in X*Y is A_k = sum_i w_i w_{i+k}."""
import json, sys
from fractions import Fraction
d = json.load(open(sys.argv[1])); K = d["K"]; w = d["weights_int"]; n = len(w)
B = 2 * max(w).bit_length() + n.bit_length() + 2
X = int.from_bytes(b"".join(x.to_bytes(B // 8 + 1, "little") for x in w), "little")  # placeholder, rebuilt below
X = 0
for i in reversed(range(n)):
    X = (X << B) | w[i]
Y = 0
for i in range(n):
    Y = (Y << B) | w[i]
Z = X * Y
mask = (1 << B) - 1
best, bk = None, None
for k in range(0, K + 1):
    A = (Z >> (B * (n - 1 + k))) & mask
    if best is None or A < best:
        best, bk = A, k
S = sum(w)
r = Fraction(K * best, S * S)
print(f"K={K} n={n} argmin_k={bk} A_min={best} S={S} score={float(r):.15f}")
