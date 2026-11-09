"""Independent check of the C1a score 2n * max_m (v*v)_m / (sum v)^2 of a non-negative integer step function.
A float64 FFT locates every lag m whose (v*v)_m is within max(1000 x error bound, 1e-10 x peak) of the largest.
Those lags are then recomputed exactly: each v_i is split into L limbs of s bits, v_i = sum_l a_{l,i} 2^{ls}, so every
partial dot product fits in int64, and the exact value is assembled in Python integers.
    python3 check1a_fft2.py file.json"""
import json, sys
from fractions import Fraction
import numpy as np

d = json.load(open(sys.argv[1])); v = [int(x) for x in d["v_int"]]; n = len(v); assert n == int(d["n"])
assert min(v) >= 0
M = max(v); L = 1
while True:
    s = -(-M.bit_length() // L)
    if 2 * s + n.bit_length() < 63: break
    L += 1
mask = (1 << s) - 1
limbs = [np.array([(x >> (l * s)) & mask for x in v], dtype=np.int64) for l in range(L)]
assert all(sum(int(limbs[l][i]) << (l * s) for l in range(L)) == v[i] for i in range(0, n, 997))

x = np.array(v, dtype=np.float64) / M
m = 1 << (2 * n - 1).bit_length()
c = np.fft.irfft(np.fft.rfft(x, m) ** 2, m)[: 2 * n - 1]
err = 2.2e-16 * np.log2(m) * float(np.dot(x, x))
tol = max(1000 * err, 1e-10 * c.max()); cand = np.nonzero(c >= c.max() - tol)[0]

def exact(k):
    lo, hi = max(0, k - n + 1), min(k, n - 1); tot = 0
    for l1 in range(L):
        p = limbs[l1][lo:hi + 1]
        for l2 in range(L):
            q = limbs[l2][k - hi:k - lo + 1][::-1]
            tot += int(np.dot(p, q)) << ((l1 + l2) * s)
    return tot

best, bk = -1, -1
for k in cand:
    val = exact(int(k))
    if val > best: best, bk = val, int(k)
lo, hi = max(0, bk - n + 1), min(bk, n - 1)
assert sum(v[i] * v[bk - i] for i in range(lo, hi + 1)) == best   # pure-Python cross-check at the peak
S = sum(v); sc = Fraction(2 * n * best, S * S)
print(f"1a fft err_bound={err:.3e} tol/max={tol / c.max():.3e} candidates={len(cand)} limbs={L}x{s}bit "
      f"peak_index={bk} peak={best}\nscore={sc.numerator}/{sc.denominator}\n~{float(sc):.15f}")
