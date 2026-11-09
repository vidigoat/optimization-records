"""Independent check: FFT autoconvolution locates every lag within the float error of the max,
then each such lag is recomputed exactly as a big-integer dot product."""
import json, sys
import numpy as np
from fractions import Fraction
d = json.load(open(sys.argv[1])); a = [int(x) for x in d["v_int"]]; n = len(a)
x = np.array(a, dtype=np.float64) / 2.0**56
m = 1 << (2 * n - 1).bit_length()
c = np.fft.irfft(np.fft.rfft(x, m) ** 2, m)[: 2 * n - 1]
# float64 FFT error is at most ~ eps * log2(m) * sum(x^2); take 1000x that, and never less than 1e-10 * max
err = 2.2e-16 * np.log2(m) * float(np.dot(x, x))
tol = max(1000 * err, 1e-10 * c.max())
cand = np.nonzero(c >= c.max() - tol)[0]
best, bk = -1, -1
for k in cand:
    k = int(k); lo, hi = max(0, k - n + 1), min(k, n - 1)
    v = sum(a[i] * a[k - i] for i in range(lo, hi + 1))
    if v > best: best, bk = v, k
S = sum(a); s = Fraction(2 * n * best, S * S)
print(f"fft err_bound={err:.3e} tol/max={tol/c.max():.3e} candidates={len(cand)} peak_index={bk} peak={best}\n~{float(s):.15f}")
