"""G(mu) = sum_{k in Z^2, k != 0} 1/(|k|^2 (1 + mu |k|^2))  =  sum_{k != 0} [1/r - 1/(r + A)],  A = 1/mu,  r = |k|^2.

Row sums (exact):  sum_y 1/(y^2 + s^2) = (pi/s) coth(pi s).
  row_0 = pi^2/3 - [(pi/a) coth(pi a) - 1/A],   a = sqrt(A)
  row_x = (pi/x) coth(pi x) - (pi/s) coth(pi s),   s = sqrt(x^2 + A),  x >= 1
G = row_0 + 2 sum_{x >= 1} row_x.
Rows x >= N: row_x = pi*phi(x) + e_x with phi(x) = 1/x - 1/sqrt(x^2+A) positive, decreasing and convex on x > 0,
so  int_N^inf phi + phi(N)/2  <=  sum_{x>=N} phi(x)  <=  int_{N-1/2}^inf phi   (trapezoid / midpoint),
int_M^inf phi = log((M + sqrt(M^2+A)) / (2M));  and 0 <= e_x <= (pi/x) * 2.0001 exp(-2 pi x) for x >= N >= 10.
"""
import math
import numpy as np
from flint import arb, ctx


def G_arb(mu, N=20000, prec=160):
    """Rigorous ball for G(mu); mu is an exact binary float (converted exactly)."""
    ctx.prec = prec
    pi = arb.pi(); A = 1 / arb(mu); a = A.sqrt()
    tot = pi * pi / 3 - (pi / a * (pi * a).coth() - 1 / A)
    acc = arb(0)
    for x in range(1, N):
        s = (arb(x * x) + A).sqrt()
        acc += pi / x * (pi * x).coth() - pi / s * (pi * s).coth()
    def tail_int(M):
        M = arb(M)
        return ((M + (M * M + A).sqrt()) / (2 * M)).log()
    phiN = 1 / arb(N) - 1 / (arb(N * N) + A).sqrt()
    lo = tail_int(N) + phiN / 2
    hi = tail_int(arb(2 * N - 1) / 2)
    ex = pi / N * arb(2.0001) * (-2 * pi * N).exp() * 2          # generous bound on sum of e_x
    tail = (lo + hi) / 2 + arb(0, ((hi - lo) / 2).upper() + ex.upper())   # ball covering [lo, hi + ex]
    return tot + 2 * (acc + pi * tail)


def G_float(mu, N=200000):
    A = 1.0 / mu; a = math.sqrt(A)
    x = np.arange(1, N, dtype=float); s = np.sqrt(x * x + A)
    rows = math.pi / x / np.tanh(math.pi * x) - math.pi / s / np.tanh(math.pi * s)
    tot = math.pi ** 2 / 3 - (math.pi / a / math.tanh(math.pi * a) - 1 / A)
    M = N - 0.5
    tail = math.log((M + math.sqrt(M * M + A)) / (2 * M))
    return tot + 2 * (np.sum(rows[::-1]) + math.pi * tail)
