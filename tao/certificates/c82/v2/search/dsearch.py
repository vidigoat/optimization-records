"""Denominator column generation by complementary-measure monic IRLS/CVP (own implementation).
rc(R) = log M(R o T) - sum_k nu_k log|R(w_k)| - lam*deg R  ~  sum_k (1/N - nu_k) log|R(w_k)| - lam*deg R.
Search R = R0 + H*U, H = X(X-1)P4, deg U < deg R0 - 6 (monic; R(0),R(1), R mod P4 fixed), IRLS on weights c_k/|R_cur(w_k)|^2.
usage: dsearch.py state_in.json state_out.json rounds starts seed"""
import numpy as np, json, sys, time, os
from multiprocessing import Pool
from fractions import Fraction
from flint import fmpz_poly, fmpz_mat, fmpq_poly, ctx
import mpmath as mp
from core import *
from wlp import WLP

N = int(os.environ.get('CG_N', 1 << 17)); t = grid(N)
X = fmpz_poly([0, 1]); P4 = fmpz_poly([1, -3, 4, -2, 1]); H = X*(X-1)*P4
DEGS = [int(x) for x in os.environ.get('DS_DEGS', '44').split(',')]
K = int(os.environ.get('DS_K', 360))
BITS = 140

def strata(nu):
    c = 1.0/N - nu; c = np.maximum(c, 0)
    cs = np.cumsum(c); tot = cs[-1]
    edges = np.searchsorted(cs, np.linspace(0, tot, K+1)[1:-1])
    groups = np.split(np.arange(N), edges)
    ts, ms = [], []
    for g in groups:
        m = c[g].sum()
        if m <= 0: continue
        ts.append((c[g] @ t[g])/m); ms.append(m)
    return np.array(ts), np.array(ms)

G = {}
def winit(nu, lam, starts_pool, nump, denp):
    ts, ms = strata(nu)
    sup = np.flatnonzero(nu > 0)
    zf = np.exp(2j*np.pi*t[sup]); G['wf'] = zf*(1-zf); G['nuf'] = nu[sup]; G['lam'] = lam
    G['ts'], G['ms'] = ts, ms
    z = np.exp(2j*np.pi*ts); G['w'] = z*(1-z)
    mp.mp.prec = 200
    wmp = []
    for tt in ts:
        zz = mp.expjpi(2*mp.mpf(tt)); wmp.append(zz*(1-zz))
    G['wmp'] = wmp
    G['pows'] = {}
    G['starts'] = [fmpz_poly(c[::-1]) for c in starts_pool]
    G['nump'] = [fmpz_poly(c[::-1]) for c in nump]
    G['denp'] = set(tuple(c) for c in denp)

def basis_evals(d):
    """integer-scaled evaluations of H*X^i, i<d-6, at the strata points: list over i of list of (re,im) mp."""
    if d not in G['pows']:
        S = mp.mpf(2)**BITS; out = []
        Hc = [int(x) for x in H.coeffs()]
        for w in G['wmp']:
            hv = mp.polyval(Hc[::-1], w); row = []; p = hv
            for i in range(d-6):
                row.append(p); p = p*w
            out.append(row)
        G['pows'][d] = out
    return G['pows'][d]

def polyval_mp(P, w):
    return mp.polyval([int(x) for x in P.coeffs()][::-1], w)

def froots(cs_desc):
    ctx.prec = 64
    out = []
    for rr, m in fmpz_poly(cs_desc[::-1]).complex_roots():
        out += [complex(rr.mid())]*m
    return np.array(out)

def rc_fast(cs_desc, r=None):
    if r is None: r = froots(cs_desc)
    lm = 0.0
    d = np.sqrt(1-4*r+0j)
    for zz in ((1+d)/2, (1-d)/2): lm += np.log(np.maximum(1.0, np.abs(zz))).sum()
    lv = np.log(np.abs(G['wf'][:, None] - r[None, :])).sum(1)
    return lm - G['nuf'] @ lv - G['lam']*(len(cs_desc)-1), lm

def lattice_step(Rc, Wt, lamreg):
    d = Rc.degree(); m = d - 6
    E = basis_evals(d); S2 = 2**BITS
    sw = [mp.sqrt(mp.mpf(float(x))) for x in Wt]
    rows = []
    M = 2**(BITS)  # anchor weight (relative to evaluation scale)
    lamI = max(1, int(lamreg * 2**BITS))
    for i in range(m):
        row = [0] + [lamI if j == i else 0 for j in range(m)]
        for k in range(len(E)):
            v = E[k][i]*sw[k]; row.append(int(mp.nint(v.real*S2))); row.append(int(mp.nint(v.imag*S2)))
        rows.append(row)
    rowT = [M] + [0]*m
    for k in range(len(E)):
        v = polyval_mp(Rc, G['wmp'][k])*sw[k]; rowT.append(int(mp.nint(v.real*S2))); rowT.append(int(mp.nint(v.imag*S2)))
    rows.append(rowT)
    L = fmpz_mat(rows).lll(delta=0.99)
    out = []
    R = [[int(L[r, c]) for c in range(m+1)] for r in range(m+1)]
    anch = [r for r in R if abs(r[0]) == M]
    short = [r for r in R if r[0] == 0][:12]
    cands = []
    for a in anch:
        s = 1 if a[0] > 0 else -1
        base = [s*x for x in a]
        cands.append(base)
        for b in short:
            cands.append([x+y for x, y in zip(base, b)]); cands.append([x-y for x, y in zip(base, b)])
    polys = []
    for c in cands:
        if any(x % lamI for x in c[1:]): continue
        U = fmpz_poly([x//lamI for x in c[1:]])
        polys.append(Rc + H*U)
    return polys

def unit_shift(R0, rng):
    a, b = int(rng.integers(-6, 7)), int(rng.integers(-6, 7))
    C = fmpq_poly(P4.coeffs()); Xq = fmpq_poly([0, 1])
    def pw(p, e):
        if e >= 0: return p**e % C
        g, s, _ = p.xgcd(C); inv = s / g
        return (inv**(-e)) % C
    target = (fmpq_poly(R0.coeffs()) * pw(Xq, a) * pw(Xq-1, b)) % C
    diff = (target - fmpq_poly(R0.coeffs())) % C
    g, s, _ = (Xq*(Xq-1)).xgcd(C); V = (diff * s / g) % C
    newp = fmpq_poly(R0.coeffs()) + Xq*(Xq-1)*V
    co = newp.numer().coeffs(); den = int(newp.denom())
    if den != 1: return R0
    return fmpz_poly([int(x) for x in co])

def job(seed):
    rng = np.random.default_rng(seed)
    d = DEGS[int(rng.integers(len(DEGS)))]
    pool = [p for p in G['starts'] if p.degree() <= d]
    R0 = fmpz_poly([1])
    for _ in range(40):
        if R0.degree() == d: break
        cand = [p for p in pool if R0.degree() + p.degree() <= d]
        if not cand: break
        R0 *= cand[int(rng.integers(len(cand)))]
    if R0.degree() != d: return []
    if rng.random() < float(os.environ.get('DS_SHIFT', 0.5)): R0 = unit_shift(R0, rng)
    if rng.random() < 0.3:
        U = fmpz_poly([int(x) for x in rng.integers(-2, 3, size=d-6)])
        R0 = R0 + H*U
    sig = float(rng.uniform(0, 1.0)) if rng.random() < 0.6 else 0.0
    pert = np.exp(sig*rng.standard_normal(len(G['ms'])))
    lamreg = 10**rng.uniform(-9, -4)
    Rc = R0; best = None; found = []
    for it in range(int(os.environ.get('DS_IT', 10))):
        cs = [int(x) for x in Rc.coeffs()][::-1]
        v = np.exp(np.log(np.abs(G['w'][:, None] - froots(cs)[None, :])).sum(1))
        Wt = G['ms']*pert/np.maximum(v, 1e-12*np.median(v))**2
        Wt = Wt/Wt.sum()
        try:
            polys = lattice_step(Rc, Wt, lamreg)
        except Exception as e:
            break
        bi = None
        for P in polys:
            c2 = [int(x) for x in P.coeffs()][::-1]
            if c2 == cs: continue
            try: r, lm = rc_fast(c2)
            except Exception: continue
            if not np.isfinite(r): continue
            if r < float(os.environ.get('DS_THR', 0.0)): found.append((r, c2))
            if bi is None or r < bi[0]: bi = (r, c2, P)
        if bi is None: break
        if best is not None and bi[0] >= best - 1e-9: break
        best = bi[0]; Rc = bi[2]
    found.sort(key=lambda x: x[0])
    return [(r, c) for r, c in found[:6]] + [('best', best, d)]

def main():
    sin, sout, rounds, starts, seed = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
    st = json.load(open(sin)); den, num = st['denominator'], st['numerator']
    L = Lib([dict(name=r['name'], coeffs=r['coeffs']) for r in den], [dict(name=s['name'], coeffs=s['coeffs']) for s in num], t)
    W = WLP(t, L.LR, L.LM, L.dR, L.LS, L.dS, thr=0.004)
    w0 = np.array(st['w'] + [0.0]*(len(num)-len(st['w']))); mu0 = np.array(st['mu'] + [0.0]*(len(den)-len(st['mu'])))
    v = W.init(mu0, w0); print('init LP', v, flush=True)
    have = {tuple(r['coeffs']) for r in den} | {tuple(s['coeffs']) for s in num}
    nump = [s['coeffs'] for s in num]
    rs = seed*1000003; nadd = 0
    for rnd in range(rounds):
        nu, price, lam = W.duals(); mu, w = W.weights()
        act = [r['coeffs'] for r, x in zip(den, mu) if x > 1e-9]
        spool = act + [r['coeffs'] for r in den if len(r['coeffs']) - 1 <= 28]
        t0 = time.time()
        near = [r['coeffs'] for r, lr_, lm_ in zip(den, W.LR, W.LM) if lm_ - nu @ lr_ - lam*(len(r['coeffs'])-1) < float(os.environ.get('DS_NEARTHR', 0.003))]
        print('   near pool', len(near), flush=True)
        with Pool(int(os.environ.get('CG_PROCS', os.cpu_count())), initializer=winit2, initargs=(nu, lam, spool, nump, [r['coeffs'] for r in den], near)) as pool:
            out = pool.map(job2, [rs+i for i in range(starts)], chunksize=1)
        rs += starts
        bests = [o[-1][1] for o in out if o and o[-1][1] is not None]
        print(f'round {rnd}: {time.time()-t0:.0f}s, starts {starts}, best-rc quantiles', np.round(np.quantile(bests, [0, .1, .5]), 6) if bests else None, flush=True)
        facs = {}
        for o in out:
            for r, cs in o[:-1]:
                for f, e in fmpz_poly(cs[::-1]).factor()[1]:
                    if f.degree() < 2: continue
                    c = [int(x) for x in f.coeffs()][::-1]
                    if c[0] < 0: c = [-x for x in c]
                    if abs(c[0]) != 1 or tuple(c) in have: continue
                    facs[tuple(c)] = c
        scored = []
        for c in facs.values():
            lr = logabs(c, t); lm = logmahlerT(c); r = lm - nu @ lr - lam*(len(c)-1)
            if r < -1e-10: scored.append((r, c, lr, lm))
        scored.sort(key=lambda x: x[0])
        print(f'   candidates {len(facs)}, negative rc {len(scored)}', [round(s[0], 6) for s in scored[:8]], flush=True)
        for r, c, lr, lm in scored[:30]:
            nu2, _, lam2 = W.duals(); r2 = lm - nu2 @ lr - lam2*(len(c)-1)
            if r2 >= -1e-12: continue
            W.add_den(lr, lm, len(c)-1); v2 = W.solve()
            if v - v2 > 1e-12:
                name = f'DS{seed}_{rnd}_{nadd}'; nadd += 1
                den.append(dict(name=name, coeffs=c)); have.add(tuple(c))
                print(f'   + {name} deg {len(c)-1} rc {r2:.4e} logM/deg {lm/(len(c)-1):.9f} gain {v-v2:.3e} -> {v2:.13f}', flush=True)
                v = v2
            else:
                W.pop_den(); W.solve()
        mu, w = W.weights()
        json.dump(dict(denominator=den, numerator=num, mu=list(mu), w=list(w), lp=v, N=N), open(sout, 'w'))


# ---------------- better CVP: LLL-reduce the H*X^i lattice, then Schnorr-Euchner enumeration around -target ----------
def enum_cvp(Rm, y, maxnodes=20000, nsol=12):
    """minimise ||y + R x||^2 over integer x, R upper triangular (n x n). returns list of (dist2, x) best nsol."""
    n = Rm.shape[0]
    diag = np.diag(Rm)
    # Babai first for radius
    x = np.zeros(n)
    for i in range(n-1, -1, -1):
        x[i] = np.round(-(y[i] + Rm[i, i+1:] @ x[i+1:]) / diag[i])
    r0 = np.sum((y + Rm @ x)**2)
    best = [(r0, x.copy())]
    radius = r0 * 1.0000001
    nodes = [0]
    xs = np.zeros(n)
    def rec(i, partial):
        if nodes[0] > maxnodes: return
        nonlocal radius
        c = -(y[i] + Rm[i, i+1:] @ xs[i+1:]) / diag[i]
        xi = np.round(c); step = 0
        # zig-zag
        k = 0
        while True:
            cand = xi + ((k+1)//2 if k % 2 else -(k//2)) * (1 if c >= xi else -1) if k > 0 else xi
            d = partial + (diag[i]*(cand - c))**2
            if d > radius:
                if k > 0: 
                    # check other side once more
                    k += 1
                    cand2 = xi + ((k+1)//2 if k % 2 else -(k//2)) * (1 if c >= xi else -1)
                    if partial + (diag[i]*(cand2 - c))**2 > radius: break
                    continue
                break
            nodes[0] += 1
            xs[i] = cand
            if i == 0:
                best.append((d, xs.copy())); best.sort(key=lambda z: z[0])
                if len(best) > nsol: best.pop()
                if len(best) == nsol: radius = best[-1][0]
            else:
                rec(i-1, d)
            if nodes[0] > maxnodes: return
            k += 1
    rec(n-1, 0.0)
    return best

def lattice_step2(Rc, Wt, maxnodes=20000, nsol=12):
    d = Rc.degree(); m = d - 6
    E = basis_evals(d); S2 = 2**BITS
    sw = [mp.sqrt(mp.mpf(float(x))) for x in Wt]
    lamI = 1
    rows = []
    for i in range(m):
        row = [lamI if j == i else 0 for j in range(m)]
        for k in range(len(E)):
            v = E[k][i]*sw[k]; row.append(int(mp.nint(v.real*S2))); row.append(int(mp.nint(v.imag*S2)))
        rows.append(row)
    L = fmpz_mat(rows).lll(delta=0.99)
    U = np.array([[int(L[r, c]) for c in range(m)] for r in range(m)], dtype=object)
    B = np.array([[float(int(L[r, c]))/float(S2) for c in range(m, L.ncols())] for r in range(m)])
    tv = []
    for k in range(len(E)):
        v = polyval_mp(Rc, G['wmp'][k])*sw[k]; tv += [float(v.real), float(v.imag)]
    tv = np.array(tv)
    Qm, Rm = np.linalg.qr(B.T)
    y = Qm.T @ tv
    sols = enum_cvp(Rm, y, maxnodes, nsol)
    out = []
    for dist, x in sols:
        xi = [int(v) for v in x]
        coeffs = [sum(xi[r]*int(U[r][c]) for r in range(m)) for c in range(m)]
        out.append((float(dist), Rc + H*fmpz_poly(coeffs)))
    return out, float(tv @ tv)

def lvq(R):
    cs = [int(x) for x in R.coeffs()][::-1]; r = froots(cs)
    return np.log(np.maximum(np.abs(G['w'][:, None] - r[None, :]), 1e-12)).sum(1), r

def job2(seed):
    rng = np.random.default_rng(seed)
    d = DEGS[int(rng.integers(len(DEGS)))]
    near = G['near']  # list of (fmpz_poly) with small rc, preferred starts
    mode = rng.random()
    if mode < float(os.environ.get('DS_PNEAR', 0.6)) and any(p.degree() == d for p in near):
        cand = [p for p in near if p.degree() == d]; R0 = cand[int(rng.integers(len(cand)))]
    else:
        pool = [p for p in G['starts'] if p.degree() <= d]
        R0 = fmpz_poly([1])
        for _ in range(40):
            if R0.degree() == d: break
            cand = [p for p in pool if R0.degree() + p.degree() <= d]
            if not cand: break
            R0 *= cand[int(rng.integers(len(cand)))]
        if R0.degree() != d: return [('best', None, d)]
    if rng.random() < float(os.environ.get('DS_SHIFT', 0.15)): R0 = unit_shift(R0, rng)
    l0, rr = lvq(R0)
    pm = rng.random()
    if pm < 0.6:
        s = 10**rng.uniform(-2.5, -0.8)
        rp = rr*np.exp(s*rng.standard_normal(len(rr))) if rng.random() < 0.5 else rr + s*np.abs(rr)*(rng.standard_normal(len(rr)) + 1j*rng.standard_normal(len(rr)))
        rp = np.concatenate([rp[np.imag(rr) >= 0], np.conj(rp[np.imag(rr) > 0])])
        logw = -2*np.log(np.maximum(np.abs(G['w'][:, None] - rp[None, :]), 1e-12)).sum(1)
    elif pm < 0.85:
        logw = -2*l0 + float(rng.uniform(0.05, 0.8))*rng.standard_normal(len(l0))
    else:
        a = rng.uniform(0.7, 1.3); logw = -2*a*l0
    Rc = R0; cur = rc_fast([int(x) for x in R0.coeffs()][::-1])[0]; found = []; best = None
    seen = {tuple(int(x) for x in R0.coeffs())}
    for it in range(int(os.environ.get('DS_IT', 8))):
        if it > 0: logw = -2*lvq(Rc)[0]
        Wt = G['ms']*np.exp(logw - logw.max()); Wt /= Wt.sum()
        try:
            sols, _ = lattice_step2(Rc, Wt, maxnodes=int(os.environ.get('DS_NODES', 4000)), nsol=int(os.environ.get('DS_NSOL', 16)))
        except Exception as e:
            break
        bi = None
        for dd, Q in sols:
            k = tuple(int(x) for x in Q.coeffs())
            if k in seen: 
                if Q == Rc and it > 0 and (bi is None or cur < bi[0]): bi = (cur, Q)
                continue
            seen.add(k)
            c2 = list(k)[::-1]
            try: r, lm = rc_fast(c2)
            except Exception: continue
            if not np.isfinite(r): continue
            if r < float(os.environ.get('DS_THR', 0.0)): found.append((r, c2))
            if bi is None or r < bi[0]: bi = (r, Q)
        if bi is None: break
        if best is None or bi[0] < best: best = bi[0]
        if it > 0 and bi[0] >= cur - 1e-9: break
        cur, Rc = bi
    found.sort(key=lambda x: x[0])
    return [(r, c) for r, c in found[:6]] + [('best', best, d)]

def winit2(nu, lam, starts_pool, nump, denp, near):
    winit(nu, lam, starts_pool, nump, denp)
    G['near'] = [fmpz_poly(c[::-1]) for c in near]

if __name__ == '__main__':
    main()
