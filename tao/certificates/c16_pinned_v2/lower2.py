"""Rigorous lower bound for C16 from the trial function u = sum_{k != 0} c_k e^{ik.x}, c_k = 1/(|k|^2 (1 + mu |k|^2)):
    L >= S0^2/(pi S1) - log delta - log(1 + log delta),  delta = S2/S1,  since ||u||_inf >= u(0) = S0.
With A = 1/mu, r = |k|^2:  S0 = G(mu),  S1 = G - A H,  S2 = A^2 H,  H = sum_{k != 0} 1/(r + A)^2.
Row sums: sum_y 1/(y^2+s^2)^2 = pi/(2 s^3) coth(pi s) + pi^2/(2 s^2) csch^2(pi s).  Rows x >= N: main term
psi(x) = pi/(2 (x^2+A)^{3/2}) is decreasing and convex for x > sqrt(A)/2, bracketed by trapezoid/midpoint;
int_M^inf psi = pi/(2A) (1 - M/sqrt(M^2+A)); the coth/csch corrections are below 4 pi exp(-2 pi N) in total.
    python lower.py mu [N]   (mu as a decimal or a hex float)"""
import sys
from flint import arb, ctx
from gsum import G_arb

def H_arb(mu, N=20000, prec=160):
    ctx.prec = prec
    pi = arb.pi(); A = 1 / arb(mu); a = A.sqrt()
    assert arb(N) - arb(1) / 2 > a / 2
    def rowsum(s): return pi / (2 * s ** 3) * (pi * s).coth() + pi * pi / (2 * s * s) / (pi * s).sinh() ** 2
    tot = rowsum(a) - 1 / (A * A)
    acc = arb(0)
    for x in range(1, N):
        acc += rowsum((arb(x * x) + A).sqrt())
    def tail_int(M):
        M = arb(M); return pi / (2 * A) * (1 - M / (M * M + A).sqrt())
    psiN = pi / (2 * (arb(N * N) + A) ** arb(1.5))
    lo = tail_int(N) + psiN / 2; hi = tail_int(arb(2 * N - 1) / 2)
    err = 4 * pi * (-2 * pi * N).exp()
    return tot + 2 * (acc + (lo + hi) / 2 + arb(0, ((hi - lo) / 2).upper() + err.upper()))

if __name__ == "__main__":
    mu = float.fromhex(sys.argv[1]) if sys.argv[1].startswith("0x") else float(sys.argv[1]); N = int(sys.argv[2]) if len(sys.argv) > 2 else 200000; ctx.prec = 192
    G = G_arb(mu, N=N, prec=192); H = H_arb(mu, N=N, prec=192); A = 1 / arb(mu); pi = arb.pi()
    S1 = G - A * H; S2 = A * A * H; d = S2 / S1
    L = G * G / (pi * S1) - d.log() - (1 + d.log()).log()
    print(f"mu = {mu!r}\ndelta = {d}\nL >= {L}\nlower endpoint {L.lower()}")
