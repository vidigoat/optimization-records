"""A lower bound for the reverse Brunn-Minkowski constant (Tao's table, entry 70a): C >= 1.0316.

Take B1 = a square of side a and B2 = a disk of radius r in the plane. For any u with |det u| = 1,
u(B1) is a parallelogram P of area a^2, and by Steiner's formula |P + B2| = a^2 + per(P) r + pi r^2.
A parallelogram with sides p, q and angle t has a^2 = p q sin t <= p q, so per(P) = 2(p + q) >=
4 sqrt(p q) >= 4a, with equality for a square. Hence min over u of |u(B1) + B2| = a^2 + 4 a r + pi r^2.
With a = sqrt(pi) r both bodies have area pi r^2, and
    |u(B1) + B2|^(1/2) / (|B1|^(1/2) + |B2|^(1/2)) >= sqrt(2 pi + 4 sqrt(pi)) / (2 sqrt(pi)) = sqrt(1/2 + 1/sqrt(pi))
for every u, so C_RBM >= sqrt(1/2 + 1/sqrt(pi)) = 1.03160... (the table lists only the trivial 1).

The script evaluates the closed form in interval arithmetic and checks the argument numerically on
random maps u (computing |u(B1) + B2| directly from the parallelogram's sides).
"""
import math
import random

from mpmath import iv, mp


def bound():
    iv.dps = 30
    return iv.sqrt(iv.mpf(1) / 2 + 1 / iv.sqrt(iv.pi))


def area_sum(u, a, r):
    """|u(square of side a) + disk of radius r|, from Steiner's formula with the exact perimeter."""
    (p, q), (s, t) = u
    e1, e2 = (p * a, s * a), (q * a, t * a)
    area = abs(e1[0] * e2[1] - e1[1] * e2[0])
    per = 2 * (math.hypot(*e1) + math.hypot(*e2))
    return area + per * r + math.pi * r * r


if __name__ == "__main__":
    b = bound()
    mp.dps = 15
    print(f"70a: C_RBM >= sqrt(1/2 + 1/sqrt(pi)) in [{mp.nstr(mp.mpf(b.a.a), 12)}, {mp.nstr(mp.mpf(b.b.b), 12)}], versus the table's trivial 1")
    r = 1.0
    a = math.sqrt(math.pi) * r
    worst = math.inf
    random.seed(1)
    for _ in range(200000):
        p, q, s = (random.uniform(-3, 3) for _ in range(3))
        if abs(p) < 1e-3:
            continue
        t = (1 + q * s) / p                       # det = p t - q s = 1
        ratio = math.sqrt(area_sum(((p, q), (s, t)), a, r)) / (a + math.sqrt(math.pi) * r)
        worst = min(worst, ratio)
    print(f"smallest ratio over 200,000 random maps with det 1: {worst:.6f} (never below the bound: {worst >= float(b.a) - 1e-12})")
