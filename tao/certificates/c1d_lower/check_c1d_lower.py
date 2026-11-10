"""Rigorous check (mpmath interval arithmetic) of the contour inequality behind  C1d >= M0.

For M = M0, eps = M - 1, the claim is: on the boundary of the rectangle [A, B] x [-C, C],
    M * |sinc(z)|  >  Bnd(|Im z|),   sinc(z) = sin(pi z)/(pi z),
    Bnd(y) = (M/(pi y)) (exp(pi y) - exp(pi y / M))   (y > 0),   Bnd(0) = M - 1,
where Bnd(y) bounds |h^(x+iy)| for every 0 <= h <= 2M on [-1/2,1/2] with integral M - 1.
Squared form checked on every box:  M^2 (sin^2(pi x) + sinh^2(pi y))  >  pi^2 (x^2 + y^2) Bnd(Ymax)^2,
using that Bnd is increasing in y.
    python3 check_c1d_lower.py [M0] [boxes_per_side]"""
import sys
from mpmath import iv, mpf
iv.prec = 120
M0 = sys.argv[1] if len(sys.argv) > 1 else "1.15"
N = int(sys.argv[2]) if len(sys.argv) > 2 else 4000
A, B, C = mpf("0.87"), mpf("1.23"), mpf("0.34")
M = iv.mpf(M0); pi = iv.pi

def bnd(Y):                      # upper bound of Bnd on [0, Y]; Y a positive mpf
    y = iv.mpf(Y)
    return (M / (pi * y)) * (iv.exp(pi * y) - iv.exp(pi * y / M))

def box_ok(x, y):                # x, y intervals
    Ymax = max(abs(mpf(y.a)), abs(mpf(y.b)))
    Ymin = mpf(0) if y.a <= 0 <= y.b else min(abs(mpf(y.a)), abs(mpf(y.b)))
    t = pi * iv.mpf(Ymin); sh = (iv.exp(t) - iv.exp(-t)) / 2      # sinh^2 is increasing in |y|
    s = iv.sin(pi * x)
    lhs = M**2 * (s**2 + sh**2)
    rhs = pi**2 * (x**2 + iv.mpf(Ymax)**2) * bnd(Ymax)**2 if Ymax > 0 else (M - 1)**2 * pi**2 * (x**2 + iv.mpf(Ymax)**2)
    return lhs.a > rhs.b, float(mpf(lhs.a) / mpf(rhs.b))

def side(p0, p1):
    worst = 1e9
    for i in range(N):
        a0 = p0[0] + (p1[0] - p0[0]) * i / N; a1 = p0[0] + (p1[0] - p0[0]) * (i + 1) / N
        b0 = p0[1] + (p1[1] - p0[1]) * i / N; b1 = p0[1] + (p1[1] - p0[1]) * (i + 1) / N
        x = iv.mpf([min(a0, a1), max(a0, a1)]); y = iv.mpf([min(b0, b1), max(b0, b1)])
        ok, r = box_ok(x, y)
        if not ok: return False, r
        worst = min(worst, r)
    return True, worst

corners = [(A, -C), (B, -C), (B, C), (A, C), (A, -C)]
allok = True
for p0, p1 in zip(corners, corners[1:]):
    ok, w = side(p0, p1); allok &= ok
    print(f"side {tuple(map(float,p0))} -> {tuple(map(float,p1))}: {'OK' if ok else 'FAIL'}  min certified ratio^2 {w:.6f}")
print(f"M0 = {M0}: contour inequality {'VERIFIED' if allok else 'NOT verified'}; the rectangle contains exactly one zero of sinc (z = 1).")
