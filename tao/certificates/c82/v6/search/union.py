"""Build union library state from several certificates. usage: union.py out.json cert1 cert2 ... (first cert gives initial weights)"""
import json, sys
from fractions import Fraction
from flint import fmpz_poly
out = sys.argv[1]; certs = [json.load(open(p)) for p in sys.argv[2:]]
den, num = [], []; seen = set()
def key(c):
    c = list(c); 
    if c[0] < 0: c = [-x for x in c]
    return tuple(c)
for ci, c in enumerate(certs):
    for r in c['denominator']:
        k = key(r['coeffs'])
        if k in seen: continue
        seen.add(k); den.append(dict(name=r['name'] + ('' if ci == 0 else f'_c{ci}'), coeffs=list(k), m0=(r['multiplicity'] if ci == 0 else 0)))
for ci, c in enumerate(certs):
    for s in c['numerator']:
        k = key(s['coeffs'])
        if k in seen: continue
        seen.add(k); num.append(dict(name=s['name'] + ('' if ci == 0 else f'_c{ci}'), coeffs=list(k), q0=(s['q_rational'] if ci == 0 else '0')))
# irreducibility / conflict check
for s in num:
    S = fmpz_poly(s['coeffs'][::-1])
    for r in den:
        assert S.gcd(fmpz_poly(r['coeffs'][::-1])).degree() == 0, (s['name'], r['name'])
D = sum((len(r['coeffs'])-1)*r['m0'] for r in den)
mu = [r['m0']/D for r in den]; w = [float(Fraction(s['q0']))/D for s in num]
json.dump(dict(denominator=[dict(name=r['name'], coeffs=r['coeffs']) for r in den], numerator=[dict(name=s['name'], coeffs=s['coeffs']) for s in num], mu=mu, w=w), open(out, 'w'))
print('den', len(den), 'num', len(num))
