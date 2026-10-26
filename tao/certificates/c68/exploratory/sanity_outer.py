# sanity: omega = z^lam * B (B annular Blaschke product with arbitrary zeros), lam chosen to make omega single-valued
import numpy as np
from sanity_enclosing import setup, Ecov, blaschke_orbit, pdist, hdist, min_enclosing, bounds
rng = np.random.default_rng(7)
c = 0.446; w, A, q = setup(c)
worst = 0; worstG = 0
for trial in range(60):
    r = rng.uniform(c*1.02, 0.98)
    nz = rng.integers(1, 3)
    phis = rng.uniform(0.03, 0.97, nz); args = rng.uniform(0, 2*np.pi, nz)
    if trial % 3 == 0:
        phis = np.r_[1 + np.log(r)/w]; args = np.r_[0.0]   # zero at r itself
    zeros = np.exp(w*(phis-1)) * np.exp(1j*args)
    def B(zeta):
        Z = Ecov(zeta, w); v = np.ones_like(Z, complex)
        for z0 in zeros: v *= blaschke_orbit(Z, Ecov(np.log(z0), w), q)
        return v
    zeta0 = np.log(r) + 0.1j
    chi = B(np.array([zeta0 + 2j*np.pi]))[0] / B(np.array([zeta0]))[0]
    gam = np.angle(chi)          # B(zeta+2 pi i) = e^{i gam} B(zeta)
    lam = (-gam/(2*np.pi)) % 1.0  # z^lam = exp(lam zeta) has character e^{2 pi i lam}
    for extra in [0, 1]:
        L = lam + extra
        th = np.linspace(0, 2*np.pi, 721)[:-1]
        zeta = np.log(r) + 1j*th
        vals = np.exp(L*zeta) * B(zeta)
        v2 = np.exp(L*(zeta0+2j*np.pi)) * B(np.array([zeta0+2j*np.pi]))[0]
        perr = abs(v2 - vals[0])
        print("max|om|", np.abs(vals).max()) if np.abs(vals).max()>=1 else None
        ratio = 0
        for i in range(0, 720, 9):
            dth = np.abs(th - th[i]); dth = np.minimum(dth, 2*np.pi - dth)
            b, GH = bounds(c, r, dth)
            d = pdist(vals, vals[i]); m = (dth > 1e-9) & (b < 1-1e-14)
            ratio = max(ratio, np.max(d[m]/b[m]))
        rho = min_enclosing(vals[::4]); _, GH = bounds(c, r, np.pi)
        Gt = (1+np.cosh(rho))/2
        worst = max(worst, ratio); worstG = max(worstG, Gt/GH)
        print('trial %d r=%.3f nz=%d lam=%.3f period_err=%.1e max d/b=%.9f G/GH=%.6f' % (trial, r, nz, L, perr, ratio, Gt/GH))
print('WORST', worst, worstG)
