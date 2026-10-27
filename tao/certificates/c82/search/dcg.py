"""Denominator column generation: R = F +- G, F a denominator factor (or product), G a product of numerator factors
with deg G < deg F (so R stays monic).  usage: dcg.py state_in.json state_out.json"""
import numpy as np, json, sys, time, itertools, os
from core import *
from wlp import WLP
from flint import fmpz_poly

N = int(os.environ.get('CG_N', 1 << 17)); t = grid(N)
st = json.load(open(sys.argv[1])); den, num = st['denominator'], st['numerator']
L = Lib([dict(name=r['name'], coeffs=r['coeffs']) for r in den], [dict(name=s['name'], coeffs=s['coeffs']) for s in num], t)
W = WLP(t, L.LR, L.LM, L.dR, L.LS, L.dS, thr=0.004)
v = W.init(np.array(st['mu']), np.array(st['w'] + [0.0]*(len(num)-len(st['w']))))
print('init', v, flush=True)
nu, price, lam = W.duals(); print('lambda', lam, flush=True)
sup = np.flatnonzero(nu > 0); ws = np.exp(2j*np.pi*t[sup]); ws = ws*(1-ws); nus = nu[sup]
P = lambda cs: fmpz_poly(cs[::-1])
Fs = [(r['name'], P(r['coeffs'])) for r in den]
# small products of denominators up to degree 32
for (a, A), (b, B) in itertools.combinations_with_replacement(Fs[:], 2):
    if A.degree() + B.degree() <= 32: Fs.append((a+'*'+b, A*B))
numP = [(s['name'], P(s['coeffs'])) for s in num]
w_, w1 = fmpz_poly([0, 1]), fmpz_poly([-1, 1])
others = [(n, p) for n, p in numP if p.degree() >= 3]
allnum = [p for _, p in numP]
have = {tuple(r['coeffs']) for r in den}
def rc_den_fast(R):
    cs = [int(x) for x in R.coeffs()][::-1]
    lm = logmahlerT(cs)
    val = np.polyval(np.array(cs, float), ws)
    return lm - nus @ np.log(np.abs(val)) - lam * R.degree(), lm
cands = {}
t0 = time.time(); ntest = 0
for fname, F in Fs:
    dF = F.degree()
    Hs = [('1', fmpz_poly([1]))] + others + [(a+'*'+b, A*B) for (a, A), (b, B) in itertools.combinations_with_replacement(others, 2) if A.degree()+B.degree() < dF]
    for hname, H in Hs:
        dh = H.degree()
        for a in range(0, dF - dh):
            for b in range(0, dF - dh - a):
                dG = dh + a + b
                if dG < dF - int(os.environ.get('DWIN', 10)) or dG >= dF: continue
                G = H * w_**a * w1**b
                for sgn in (1, -1):
                    R = F + sgn*G
                    ntest += 1
                    for f, e in R.factor()[1]:
                        if f.degree() < 2: continue
                        cs = [int(x) for x in f.coeffs()][::-1]
                        if cs[0] < 0: cs = [-x for x in cs]
                        if abs(cs[0]) != 1 or tuple(cs) in have or tuple(cs) in cands: continue
                        r, lm = rc_den_fast(f)
                        cands[tuple(cs)] = (r, lm, f'{fname}{"+" if sgn>0 else "-"}{hname}*w^{a}*(w-1)^{b}')
print('tested', ntest, 'distinct factors', len(cands), 'time %.0f' % (time.time()-t0), flush=True)
neg = sorted([(x[0], k, x[1], x[2]) for k, x in cands.items() if x[0] < 0])
print('negative rc', len(neg), flush=True)
for r, k, lm, src in neg[:15]: print('  rc %.3e  logM/deg %.5f deg %d  %s' % (r, lm/(len(k)-1), len(k)-1, src), flush=True)
added = []
for r, k, lm, src in neg[:200]:
    R = P(list(k))
    if any(R.gcd(p).degree() > 0 for p in allnum): continue
    cs = list(k); lr = logabs(cs, t)
    nu2, _, lam2 = W.duals(); r2 = lm - nu2 @ lr - lam2*(len(cs)-1)
    if r2 >= -1e-12: continue
    W.add_den(lr, logmahlerT(cs), len(cs)-1); v2 = W.solve()
    if v - v2 > 1e-13:
        name = 'D%d' % len(added); added.append(name)
        den.append(dict(name=name, coeffs=cs, source=src)); print(f'   + {name} deg {len(cs)-1} rc {r2:.3e} gain {v-v2:.3e} -> {v2:.13f}  ({src})', flush=True)
        v = v2
    else:
        W.pop_den(); W.solve()
mu, w = W.weights()
json.dump(dict(denominator=den, numerator=num, mu=list(mu), w=list(w), lp=v, N=N), open(sys.argv[2], 'w'))
print('final', v)
