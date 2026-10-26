"""Exploratory: merge split atom pairs and re-optimise node positions (not part of proof)."""
import numpy as np, json, argparse
from fractions import Fraction as F
import build_cert as bc

def merge(nodes, wts, gap):
    # cluster consecutive nodes whose log(1-t) differ by < gap
    u = np.log(1 - nodes); out_n = []; out_w = []
    i = 0
    while i < len(nodes):
        j = i
        while j + 1 < len(nodes) and u[j + 1] - u[j] > -gap and abs(u[j+1]-u[i]) < gap: j += 1
        W = wts[i:j + 1].sum()
        if i == 0 or j == len(nodes) - 1:      # keep endpoints fixed (c^2 and R)
            for k in range(i, j + 1):
                if wts[k] > 0: out_n.append(nodes[k]); out_w.append(wts[k])
        else:
            uu = (u[i:j + 1] * wts[i:j + 1]).sum() / W
            out_n.append(1 - np.exp(uu)); out_w.append(W)
        i = j + 1
    return np.array(out_n), np.array(out_w)

ap = argparse.ArgumentParser()
ap.add_argument('--c', required=True); ap.add_argument('--R', required=True)
ap.add_argument('--kmax', type=int, default=6000); ap.add_argument('--out', required=True)
a = ap.parse_args()
cF, RF = F(a.c), F(a.R); c, R = float(cF), float(RF)
ks = np.arange(1, a.kmax + 1); pts = bc.partition(R)
s = c * c
u = np.linspace(np.log(1 - s), np.log(1 - R), 1200)
nodes = 1 - np.exp(u); nodes[0] = s; nodes[-1] = R
res = bc.lp(c, R, nodes, ks, pts, wlast_min=R / a.kmax)
wts = res.x[:-1]; keep = wts > 1e-14; keep[0] = keep[-1] = True
nodes, wts = nodes[keep], wts[keep]
print('grid tau', res.x[-1], len(nodes))
nodes, wts = merge(nodes, wts, 0.05)
for h in [2e-2, 1e-2, 5e-3, 2e-3, 1e-3, 5e-4, 2e-4, 1e-4, 5e-5, 2e-5, 1e-5]:
    for rep in range(3):
        cand = [nodes[0]]
        for t in nodes[1:-1]:
            ul = np.log(1 - t)
            cand += [1 - np.exp(ul + h * m) for m in (-2, -1, 0, 1, 2)]
        cand += [nodes[-1]]
        cand = np.array(sorted(set(x for x in cand if x <= R)))
        r = bc.lp(c, R, cand, ks, pts, wlast_min=R / a.kmax)
        w = r.x[:-1]; kk = w > 1e-15; kk[0] = kk[-1] = True
        nodes, wts = merge(cand[kk], w[kk], 4.5 * h)
    r = bc.lp(c, R, nodes, ks, pts, wlast_min=R / a.kmax)
    print('h=%g atoms=%d tau(merged)=%.4e' % (h, len(nodes), r.x[-1]))
wts = r.x[:-1]
print('FINAL nodes', repr(nodes)); print('FINAL w', repr(wts), r.x[-1])
D = 10**18
tq = [cF * cF] + [F(int(round(x * D)), D) for x in nodes[1:-1]] + [RF]
wq = [F(int(round(x * D)), D) for x in wts]
wq[-1] = max(wq[-1], F(int(round(R / a.kmax * D)) + 1, D))
j = int(np.argmax(wts)); wq[j] += RF - sum(wq)
assert sum(wq) == RF and all(x > 0 for x in wq)
Am, Bm = bc.margins(c, R, [float(x) for x in tq], [float(x) for x in wq], ks, pts)
print('float margins', Am.min(), Bm.min())
npts = len(pts) - 1
tp = [1 - (1 - RF) * F(4, 5)**i for i in range(npts)] + [F(1)]
json.dump(dict(model='flux_halfdiameter', c=str(cF), R=str(RF), nodes=[str(x) for x in tq],
               weights=[str(x) for x in wq], tail_partition=[str(x) for x in tp]), open(a.out, 'w'), indent=1)
