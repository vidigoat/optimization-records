"""Numerical cross-checks for the bound C41b <= 2*sqrt(2) - 1 (the proof is in README.md; this script is a sanity check).

1. The (s,t)-coordinate description of H ∩ (left corridor at 45°) ∩ (right corridor at 45°) agrees with direct
   polygon geometry (shapely) on random positions.
2. Random search over the three free parameters never exceeds 2*sqrt(2) - 1, and gets close to it.
3. The two one-variable inequalities used in the proof hold (checked on a fine grid; they are also proved in README).
    python3 check_sofa.py          (needs numpy and shapely)"""
import numpy as np
from shapely.geometry import Polygon, box
from shapely import affinity

r = np.sqrt(2); BIG = 200.0; Lr = 40.0
W1 = Polygon([(-Lr, 0), (1, 0), (1, Lr), (0, Lr), (0, 1), (-Lr, 1)])     # arm from the left, turning up (left turn)
W2 = Polygon([(0, -Lr), (1, -Lr), (1, 0), (Lr, 0), (Lr, 1), (0, 1)])     # arm from below, turning right (right turn)
H = box(-Lr, 0, Lr, 1)

def corridor(W, deg, x):           # the corridor as seen from the sofa, when the sofa's orientation is `deg`
    return affinity.translate(affinity.rotate(W, -deg, origin=(0, 0)), x[0], x[1])

def F_st(a, b, c, d):              # area, in (s,t) = (x-y, x+y) coordinates, of D ∩ L1 ∩ L2 (README notation)
    D = Polygon([(-BIG, -BIG), (BIG, BIG), (BIG, BIG + 2), (-BIG, -BIG + 2)])
    L1 = box(a, -BIG, a + r, b).union(box(a, b - r, BIG, b))
    L2 = box(-BIG, c, d, c + r).union(box(d - r, c, d, BIG))
    return D.intersection(L1).intersection(L2).area

rng = np.random.default_rng(0)
err = 0.0
for _ in range(500):
    x1, x2 = rng.uniform(-2.5, 2.5, 2), rng.uniform(-2.5, 2.5, 2)
    direct = H.intersection(corridor(W1, 45, x1)).intersection(corridor(W2, 45, x2)).area
    s1, t1, s2, t2 = x1[0] - x1[1], x1[0] + x1[1], x2[0] - x2[1], x2[0] + x2[1]
    err = max(err, abs(F_st(s2, r + t2, t1, r + s1) / 2 - direct))
print(f"1. max |(s,t) formula - direct geometry| over 500 random positions: {err:.2e}")

best = 0.0
for _ in range(30000):
    best = max(best, F_st(0.0, rng.uniform(-1, 7), rng.uniform(-4, 4), rng.uniform(0, 6)) / 2)
print(f"2. largest area found by random search: {best:.10f}   bound 2*sqrt(2)-1 = {2 * r - 1:.10f}")

w = np.linspace(0, r, 200001)
g = np.where(w <= 2 - r, w * r, w * r - (w + r - 2) ** 2 / 4)
print(f"3. min over w in [0, sqrt 2] of (sqrt2 - 1 + w) - g(w) = {np.min(r - 1 + w - g):.3e}  (must be >= 0; 0 at w = sqrt 2)")
