# numeric sanity check of the flux lemma constants against the explicit degree-2 extremal
import mpmath as mp
mp.mp.dps = 50
c = mp.mpf('0.44'); w = -mp.log(c); A = mp.pi**2 / w
def E(zeta): return mp.exp(1j*mp.pi*(zeta + w)/w)
def G(z, a, K=40):
    zz = E(mp.log(z)); s = 0
    for k in range(-K, K+1):
        aa = E(mp.log(a) + 2j*mp.pi*k)
        s += -mp.log(abs((zz-aa)/(zz-mp.conj(aa))))
    return s
for r in [mp.mpf(x) for x in ['0.45','0.5','0.6','0.66','0.75','0.9','0.97']]:
    phi = 1 + mp.log(r)/w; s2 = mp.sin(mp.pi*phi)**2
    S = mp.sinh(A/2)**2
    KP = S/(S+s2)
    g1 = G(-r, r); g2 = G(-r, c/r)
    bound_term = s2/mp.cosh(A/2)**2
    # also check other second zeros on the same radius at other angles (should be larger)
    g2b = min(G(-r, (c/r)*mp.expj(th)) for th in mp.linspace(-0.5,0.5,11))
    print(float(r), 'G(-r,r)=%.6e  -logKP=%.6e  G(-r,c/r)=%.6e  bound=%.6e ratio=%.8f min_nearby=%.6e' % (g1, -mp.log(KP), g2, bound_term, g2/bound_term, g2b))
