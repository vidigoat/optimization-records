"""Fine-grid LP + exact rational rounding -> certificate JSON (Pinho's format).
usage: finalize.py state.json out_cert.json [Nfine_log2=19] [Neval_log2=20]"""
import numpy as np, json, sys
from fractions import Fraction
from core import *
from wlp import WLP

st = json.load(open(sys.argv[1])); out = sys.argv[2]
nf = int(sys.argv[3]) if len(sys.argv) > 3 else 19
ne = int(sys.argv[4]) if len(sys.argv) > 4 else 20
den, num = st['denominator'], st['numerator']
t = grid(1 << nf)
L = Lib([dict(name=r['name'], coeffs=r['coeffs']) for r in den], [dict(name=s['name'], coeffs=s['coeffs']) for s in num], t)
mu0, w0 = np.array(st['mu']), np.array(st['w'])
W = WLP(t, L.LR, L.LM, L.dR, L.LS, L.dS, thr=0.003, budget=1-1e-6)
v = W.init(mu0, w0); mu, w = W.weights()
print('fine LP', nf, v, flush=True)
Dscale = 10**12
m = [max(1, int(round(x * Dscale))) for x in mu]
D = sum(mi * (len(r['coeffs'])-1) for mi, r in zip(m, den))
numout = []
for s, x in zip(num, w):
    if x * D < 1e-3: continue
    q = Fraction(int(x * D * (1 - 1e-9) * 2**20), 2**20)
    if q > 0: numout.append(dict(name=s['name'], coeffs=s['coeffs'], q_rational=str(q)))
degA = sum((len(s['coeffs'])-1) * Fraction(s['q_rational']) for s in numout)
assert degA < D, (degA, D)
cert = dict(description="Doche limit-point construction (C82, Zhang-Zagier height). Denominator factors and the first 31 numerator factors from R. Pinho's certificate (github.com/rpinho/c82-height-certificate); additional numerator factors found by IRLS dual-weighted lattice search with LLL; weights jointly re-optimised by LP and rounded exactly.",
            target=None, grid=1 << nf, numeric_unverified=None,
            denominator=[dict(name=r['name'], coeffs=r['coeffs'], multiplicity=mi) for r, mi in zip(den, m)],
            numerator=numout,
            credits=["Doche: limit-point construction and original polynomial library", "Grinsztajn: perturbative denominator factors", "Epi/BakeLab: Ed5 denominator and unequal multiplicities", "R. Pinho: 25 numerator factors and certificate format"])
# numeric evaluation on a finer grid
te = grid(1 << ne)
Dn = D
LQ = sum(r['multiplicity'] * logabs(r['coeffs'], te) for r in cert['denominator']) / Dn
LA = sum(float(Fraction(s['q_rational'])) * logabs(s['coeffs'], te) for s in numout) / Dn
Mh = sum(r['multiplicity'] * logmahlerT(r['coeffs']) for r in cert['denominator']) / Dn
val = Mh + np.maximum(LA - LQ, 0).mean()
print('rounded numeric value (grid 2^%d): %.13f' % (ne, val), 'numerators', len(numout), 'D', D, 'degA/D', float(degA / D), flush=True)
cert['numeric_unverified'] = float(val)
cert['target'] = '%.10f' % (np.ceil(val * 1e10 + 1) / 1e10)
json.dump(cert, open(out, 'w'), indent=1)
