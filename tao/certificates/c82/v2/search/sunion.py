"""union of search states: base state keeps its mu/w; other states contribute new factors with weight 0.
usage: sunion.py out.json base.json other1.json ..."""
import json, sys
from flint import fmpz_poly
out, base = sys.argv[1], json.load(open(sys.argv[2]))
den, num = base['denominator'], base['numerator']; mu, w = list(base['mu']), list(base['w'])
w += [0.0]*(len(num)-len(w)); mu += [0.0]*(len(den)-len(mu))
have = {tuple(r['coeffs']) for r in den} | {tuple(s['coeffs']) for s in num}
for p in sys.argv[3:]:
    o = json.load(open(p))
    for r in o['denominator']:
        if tuple(r['coeffs']) not in have:
            R = fmpz_poly(r['coeffs'][::-1])
            if any(R.gcd(fmpz_poly(s['coeffs'][::-1])).degree() > 0 for s in num): print('skip den', r['name']); continue
            den.append(dict(name=r['name'], coeffs=r['coeffs'])); mu.append(0.0); have.add(tuple(r['coeffs'])); print('add den', r['name'], len(r['coeffs'])-1)
    for s in o['numerator']:
        if tuple(s['coeffs']) not in have:
            S = fmpz_poly(s['coeffs'][::-1])
            if any(S.gcd(fmpz_poly(r['coeffs'][::-1])).degree() > 0 for r in den): print('skip num', s['name']); continue
            num.append(dict(name=s['name'], coeffs=s['coeffs'])); w.append(0.0); have.add(tuple(s['coeffs']))
json.dump(dict(denominator=den, numerator=num, mu=mu, w=w, lp=base.get('lp')), open(out, 'w'))
print('den', len(den), 'num', len(num))
