"""prune a state: solve LP, drop columns with zero weight and reduced cost > thr. usage: prune.py in out [thr]"""
import numpy as np, json, sys
from core import *
from wlp import WLP
st = json.load(open(sys.argv[1])); thr = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-3; thrd = float(sys.argv[4]) if len(sys.argv) > 4 else thr
t = grid(1 << 17); den, num = st['denominator'], st['numerator']
L = Lib([dict(name=r['name'], coeffs=r['coeffs']) for r in den], [dict(name=s['name'], coeffs=s['coeffs']) for s in num], t)
W = WLP(t, L.LR, L.LM, L.dR, L.LS, L.dS, thr=0.004)
mu0 = np.array(st['mu'] + [0.0]*(len(den)-len(st['mu']))); w0 = np.array(st['w'] + [0.0]*(len(num)-len(st['w'])))
v = W.init(mu0, w0); mu, w = W.weights(); nu, price, lam = W.duals()
rcd = L.LM - L.LR @ nu - lam*L.dR; rcn = L.LS @ nu + price*L.dS
kd = [i for i in range(len(den)) if mu[i] > 1e-12 or rcd[i] < thrd]
kn = [i for i in range(len(num)) if w[i] > 1e-12 or rcn[i] < thr]
print('rc den sorted', np.round(np.sort(rcd)[:30], 4)); print('LP', v, 'keep den', len(kd), '/', len(den), 'num', len(kn), '/', len(num))
json.dump(dict(denominator=[den[i] for i in kd], numerator=[num[i] for i in kn], mu=[float(mu[i]) for i in kd], w=[float(w[i]) for i in kn], lp=v), open(sys.argv[2], 'w'))
