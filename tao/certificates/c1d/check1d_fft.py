"""Independent check of the C1d score 2n * max_m |(v*v)_m| / (sum v)^2 of a signed integer step function.
A float64 FFT locates every lag m whose |(v*v)_m| is within max(100 x error bound, 1e-11 x peak) of the largest.
Those lags are then recomputed exactly: each v_i is split as v_i = a_i 2^s + b_i (0 <= b_i < 2^s), so every partial
dot product fits in int64 and the exact value is assembled in Python integers.
    python3 check1d_fft.py file.json"""
import json, sys
from fractions import Fraction
import numpy as np

d = json.load(open(sys.argv[1])); v = [int(x) for x in d["v_int"]]; n = len(v); assert n == int(d["n"])
M = max(abs(x) for x in v); s = (M.bit_length() + 1) // 2
assert 2 * s + n.bit_length() < 63 and 2 * (M.bit_length() - s + 1) + n.bit_length() < 63
A = np.array([x >> s for x in v], dtype=np.int64); B = np.array([x & ((1 << s) - 1) for x in v], dtype=np.int64)
assert all(int(a) * (1 << s) + int(b) == x for a, b, x in zip(A, B, v))

x = np.array(v, dtype=np.float64) / M
m = 1 << (2 * n - 1).bit_length()
c = np.abs(np.fft.irfft(np.fft.rfft(x, m) ** 2, m)[: 2 * n - 1])
err = 2.2e-16 * np.log2(m) * float(np.dot(x, x))
tol = max(100 * err, 1e-11 * c.max()); cand = np.nonzero(c >= c.max() - tol)[0]

def exact(k):
    lo, hi = max(0, k - n + 1), min(k, n - 1)
    a1, b1 = A[lo:hi + 1], B[lo:hi + 1]; a2, b2 = A[k - hi:k - lo + 1][::-1], B[k - hi:k - lo + 1][::-1]
    aa, ab, ba, bb = (int(np.dot(p, q)) for p, q in ((a1, a2), (a1, b2), (b1, a2), (b1, b2)))
    return (aa << (2 * s)) + ((ab + ba) << s) + bb

best, bk = -1, -1
for k in cand:
    val = abs(exact(int(k)))
    if val > best: best, bk = val, int(k)
# spot-check the int64 splitting against a pure-Python dot product at the peak
lo, hi = max(0, bk - n + 1), min(bk, n - 1)
assert abs(sum(v[i] * v[bk - i] for i in range(lo, hi + 1))) == best
S = sum(v); sc = Fraction(2 * n * best, S * S)
print(f"1d fft err_bound={err:.3e} tol/max={tol / c.max():.3e} candidates={len(cand)} peak_index={bk} |peak|={best}\n"
      f"score={sc.numerator}/{sc.denominator}\n~{float(sc):.15f}")
