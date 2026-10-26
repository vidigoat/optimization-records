"""Exploratory certificate builder (floating point; NOT part of the proof).
Model 'H': g(t) = 2 sqrt(2 s)/(sqrt(2 s) + sqrt(2 s + Psi)),  s = sin^2(pi log t/(2w)),
Psi = (1+q)^2/(8q), q = exp(-pi^2/w), w = -log c.
Tail [R,1]: piecewise-linear lower bound l_i (1-t) on a geometric partition.
Outputs a JSON certificate with exact rationals."""
import numpy as np, json, argparse
from fractions import Fraction as F
from scipy.optimize import linprog

def model(c):
    w = -np.log(c); A = np.pi**2 / w; q = np.exp(-A); Psi = (1 + q)**2 / (8 * q)
    def g(t):
        s = np.sin(np.pi * np.log(t) / (2 * w))**2
        return 2 * np.sqrt(2 * s) / (np.sqrt(2 * s) + np.sqrt(2 * s + Psi))
    return w, A, Psi, g

def piece_slopes(c, pts):
    """l_i for piece [pts[i], pts[i+1]] (floats), same formula as the checker."""
    w, A, Psi, g = model(c)
    out = []
    for a, b in zip(pts[:-1], pts[1:]):
        dhi, dlo = 1 - a, 1 - b
        U = -np.pi * np.log(a) / (2 * w)
        bb = (np.pi / (2 * w)) * (1 + dlo / 2) * (1 - U * U / 6)
        l = 2 * np.sqrt(2) * bb / (np.sqrt(2) * bb * dhi + np.sqrt(2 * bb * bb * dhi * dhi + Psi))
        out.append(l)
    return np.array(out)

def tail_lb(c, pts, ks):
    l = piece_slopes(c, pts)
    pts = np.asarray(pts, float)
    res = np.zeros(len(ks))
    for i, k in enumerate(ks):
        J = pts**(k + 1) / (k + 1) - pts**(k + 2) / (k + 2)
        res[i] = np.sum(l * (J[1:] - J[:-1]))
    return res

def partition(R, ratio=0.8, dmin=1e-8):
    d = [1 - R]
    while d[-1] * ratio > dmin:
        d.append(d[-1] * ratio)
    pts = [1 - x for x in d] + [1.0]
    return pts

def lp(c, R, nodes, ks, pts, fixed_w=None, wlast_min=0.0):
    w, A, Psi, g = model(c)
    nodes = np.asarray(nodes, float)
    beta = 1 - g(nodes)
    beta[np.abs(nodes - c * c) < 1e-15] = 1.0
    I = tail_lb(c, pts, ks)
    lsk = (ks + 1) * np.log(R) - np.log(ks + 1)
    P = np.exp(np.outer(ks, np.log(nodes)) - lsk[:, None])
    n = len(nodes)
    A1 = np.hstack([-P, np.ones((len(ks), 1))]); b1 = -np.ones(len(ks))
    A2 = np.hstack([P * beta[None, :], np.ones((len(ks), 1))]); b2 = I / np.exp(lsk) + 1
    res = linprog(np.r_[np.zeros(n), -1], A_ub=np.vstack([A1, A2]), b_ub=np.r_[b1, b2],
                  A_eq=np.r_[np.ones(n), 0][None, :], b_eq=[R],
                  bounds=[(0, None)] * (n - 1) + [(wlast_min, None)] + [(None, 1)], method='highs')
    return res

def margins(c, R, nodes, wts, ks, pts):
    w, A, Psi, g = model(c)
    nodes = np.asarray(nodes, float); wts = np.asarray(wts, float)
    beta = 1 - g(nodes); beta[np.abs(nodes - c * c) < 1e-15] = 1.0
    I = tail_lb(c, pts, ks)
    sk = np.exp((ks + 1) * np.log(R) - np.log(ks + 1))
    Pw = np.exp(np.outer(ks, np.log(nodes)))
    Am = (Pw @ wts - sk) / sk
    Bm = (I + sk - Pw @ (wts * beta)) / sk
    return Am, Bm

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--c', required=True)          # rational string e.g. 223/500
    ap.add_argument('--R', required=True)          # rational string e.g. 99/100
    ap.add_argument('--ngrid', type=int, default=1200)
    ap.add_argument('--kmax', type=int, default=6000)
    ap.add_argument('--out', required=True)
    ap.add_argument('--wtol', type=float, default=1e-13)
    a = ap.parse_args()
    cF, RF = F(a.c), F(a.R)
    c, R = float(cF), float(RF)
    s = c * c
    u = np.linspace(np.log(1 - s), np.log(1 - R), a.ngrid)
    nodes = 1 - np.exp(u); nodes[0] = s; nodes[-1] = R
    ks = np.arange(1, a.kmax + 1)
    pts = partition(R)
    res = lp(c, R, nodes, ks, pts, wlast_min=R / a.kmax)
    print('LP status', res.status, 'tau', res.x[-1])
    wts = res.x[:-1]
    idx = np.where(wts > a.wtol)[0]
    if 0 not in idx: idx = np.r_[0, idx]
    if len(nodes) - 1 not in idx: idx = np.r_[idx, len(nodes) - 1]
    nd = nodes[idx]
    # re-solve LP on the support only (cleans tiny weights)
    res2 = lp(c, R, nd, ks, pts, wlast_min=R / a.kmax)
    print('support LP', res2.status, 'tau', res2.x[-1], 'atoms', len(nd))
    wt = res2.x[:-1]
    keep = wt > 1e-16
    keep[0] = True; keep[-1] = True
    nd = nd[keep]; wt = wt[keep]
    # rationalise
    D = 10**18
    tq = [cF * cF] + [F(int(round(x * D)), D) for x in nd[1:-1]] + [RF]
    wq = [F(int(round(x * D)), D) for x in wt]
    wq[-1] = max(wq[-1], F(int(round(R / a.kmax * D)) + 1, D))
    j = int(np.argmax(wt)); wq[j] += RF - sum(wq)
    assert sum(wq) == RF and all(x > 0 for x in wq)
    Am, Bm = margins(c, R, [float(x) for x in tq], [float(x) for x in wq], ks, pts)
    print('float margins: minA %.3e at k=%d, minB %.3e at k=%d' % (Am.min(), ks[Am.argmin()], Bm.min(), ks[Bm.argmin()]))
    # exact partition: 1-(1-R)(4/5)^i, i=0..n-1, then 1
    npts = len(pts) - 1
    tp = [1 - (1 - RF) * F(4, 5)**i for i in range(npts)] + [F(1)]
    cert = dict(model='flux_halfdiameter', c=str(cF), R=str(RF),
                nodes=[str(x) for x in tq], weights=[str(x) for x in wq],
                tail_partition=[str(x) for x in tp])
    json.dump(cert, open(a.out, 'w'), indent=1)
    print('wrote', a.out)
