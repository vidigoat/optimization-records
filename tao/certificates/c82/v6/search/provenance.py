"""Classify the factors of a C82 certificate by origin, by matching full coefficient arrays.
usage: provenance.py cert.json pr330_cert.json our_cert1.json [our_cert2.json ...]"""
import json, sys
key = lambda cs: tuple(cs) if cs[0] > 0 else tuple(-x for x in cs)
c = json.load(open(sys.argv[1])); p = json.load(open(sys.argv[2]))
pin = {key(r['coeffs']) for r in p['denominator'] + p['numerator']}
ours = {}
for f in sys.argv[3:]:
    d = json.load(open(f))
    for r in d['denominator'] + d['numerator']: ours.setdefault(key(r['coeffs']), f)
for part in ('denominator', 'numerator'):
    src = {}
    for r in c[part]:
        k = key(r['coeffs']); n = r['name']
        if n in ('P1', 'P2', 'P3', 'P4', 'P6', 'P7'): s = 'Doche'
        elif k in pin: s = 'PR #330'
        elif k in ours: s = 'ours, ' + ours[k]
        else: s = 'new in this certificate'
        src.setdefault(s, []).append(n)
    print(part, len(c[part]))
    for s, v in src.items(): print('  ', s, len(v), v)
print('non-monic numerators', sum(abs(s['coeffs'][0]) != 1 for s in c['numerator']))
