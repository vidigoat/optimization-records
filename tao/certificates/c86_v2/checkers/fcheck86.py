"""Float-level sanity check of a C86 certificate (not a proof): bound, degree condition, integrals, f on a fine grid.
    python fcheck86.py cert.json [points]"""
import sys, json
from fractions import Fraction as F
import numpy as np
sys.path.insert(0, '.')
from pl import hat_potentials, poly_logabs, poly_integral_rows

c = json.load(open(sys.argv[1])); npts = int(sys.argv[2]) if len(sys.argv) > 2 else 400001
e = np.array([float(F(v)) for v in c['edges']]); h = np.array([float(F(v)) for v in c['node_densities']])
b = np.array([float(F(v)) for v in c['b']]); A = c['A']
assert h[0] == 0 and h[-1] == 0 and np.all(h >= 0) and np.all(np.diff(e) > 0) and e[0] > 0
hi = h[1:-1]
mass = np.sum(np.diff(e) * (h[1:] + h[:-1]) / 2)
print('bound %.13f  mass %.12f  sum b deg/m %.6f' % (float(F(c['bound'])), mass, sum(bj * (len(a) - 1) for a, bj in zip(A, b)) / mass))
I = np.array([poly_integral_rows(a, e) @ hi for a in A]); print('min integral %.3e' % I.min())
sup = np.nonzero(h > 0)[0]; lo, hiE = e[sup[0] - 1], e[sup[-1] + 1]
x = np.linspace(lo, hiE, npts); f = np.zeros(len(x))
for i in range(0, len(x), 3000): f[i:i + 3000] = hat_potentials(e, x[i:i + 3000]) @ hi
f -= sum(bj * poly_logabs(a, x) for a, bj in zip(A, b))
print('support [%.5f, %.5f], min f on %d points: %.3e at %.6f' % (lo, hiE, npts, f.min(), x[f.argmin()]))
