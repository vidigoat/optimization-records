"""Column generation: LP (warm) + IRLS lattice search for new numerator factors.
usage: cg.py state_in.json state_out.json rounds starts_per_round seed"""
import numpy as np, json, sys, time, os
from multiprocessing import Pool
from core import *
from wlp import WLP
from irls import Support, lattice
from flint import fmpz_poly

N = int(os.environ.get('CG_N', 1 << 17))
t = grid(N)

def irr_factors(cs):
    out = []
    for f, e in fmpz_poly(cs[::-1]).factor()[1]:
        c = [int(x) for x in f.coeffs()][::-1]
        if c[0] < 0: c = [-x for x in c]
        if len(c) >= 2: out.append(c)
    return out

G = {}
def winit(nu, den_coeffs, lib_coeffs, wpos):
    G['sp'] = Support(t, nu, sub=int(os.environ.get('CG_SUB', 3)))
    G['dens'] = [fmpz_poly(c[::-1]) for c in den_coeffs]
    G['lib'] = [fmpz_poly(c[::-1]) for c in lib_coeffs]
    G['wpos'] = wpos
    from gen2 import clusters
    G['clusters'] = clusters(t, nu, None)

def ok(P):
    return P.degree() >= 1 and all(P.gcd(x).degree() == 0 for x in G['dens'])

def irls_job(args):
    seed, = args
    rng = np.random.default_rng(seed)
    sp = G['sp']; lib = G['lib']; wpos = G['wpos']
    p = wpos / wpos.sum(); p = 0.5*p + 0.5/len(lib)
    mode = rng.random()
    if mode < 0.5:
        # real start polynomial with zeros at random dual atoms (clusters weighted by mass)
        cl = G['clusters']
        K = int(rng.integers(3, int(os.environ.get('CG_KMAX', 26))))
        pc = cl[:, 1] ** float(rng.uniform(0.3, 1.0)); pc /= pc.sum()
        pick = rng.choice(len(cl), size=min(K, len(cl)), replace=False, p=pc)
        z = np.exp(2j*np.pi*cl[pick, 0]); a_ = z*(1-z)
        a_ = a_ + (rng.standard_normal(len(a_)) + 1j*rng.standard_normal(len(a_))) * float(10**rng.uniform(-4, -1.5))
        rts = np.concatenate([a_, a_.conj()])
        c0 = np.real(np.poly(rts))
        if rng.random() < 0.5:
            f = lib[int(rng.choice(len(lib), p=p))]
            c0 = np.convolve(c0, np.array([int(x) for x in f.coeffs()][::-1], float))
        cs0f = c0
        cs0 = None
    else:
        D0 = int(rng.integers(10, int(os.environ.get('CG_DMAX', 62))))
        S0 = fmpz_poly([1])
        if rng.random() < 0.3:
            S0 *= G['dens'][int(rng.integers(len(G['dens'])))]
        for _ in range(12):
            f = lib[int(rng.choice(len(lib), p=p))]
            if S0.degree() + f.degree() <= D0: S0 *= f
        if S0.degree() < 3: return []
        cs0 = [int(x) for x in S0.coeffs()][::-1]
        cs0f = np.array(cs0, float)
    sig = float(rng.uniform(0, float(os.environ.get('CG_SIG', 1.0)))) if rng.random() < 0.6 else 0.0
    pert = np.exp(sig * rng.standard_normal(len(sp.nu)))
    extra = int(rng.integers(int(os.environ.get('CG_EXMIN', 0)), 9)); lam_rel = 10**rng.uniform(-6, -1); nk = int(rng.integers(6, 16))
    best_ok = np.inf
    found = []
    for it in range(20):
        d = len(cs0f) - 1 + extra
        if d > int(os.environ.get('CG_DCAP', 80)): break
        v = np.abs(np.polyval(cs0f, sp.w))
        fl = np.median(v) * 1e-6
        W = sp.nu * pert / np.maximum(v, fl)**2
        lam = lam_rel * np.sqrt(sp.nu.sum()) / np.linalg.norm(cs0f) * (1 if cs0 is not None else float(10**rng.uniform(0, 4)))
        try:
            B = lattice(sp, d, W, lam, bits=160)
        except Exception:
            break
        cands = []
        for i, b in enumerate(B[:nk]):
            cands.append(b)
            for b2 in B[:i]:
                cands.append([x+y for x, y in zip(b, b2)]); cands.append([x-y for x, y in zip(b, b2)])
        best = None
        for q in cands:
            P = fmpz_poly(q)
            if P.degree() < 1: continue
            cs = [int(x) for x in P.coeffs()][::-1]
            if cs[0] < 0: cs = [-x for x in cs]
            if cs == cs0 or not ok(P): continue
            r = sp.rc_fast(cs)
            if r < float(G.get("thr", -1e-10)): found.append((r, cs))
            if best is None or r < best[0]: best = (r, cs)
        if best is None or best[0] >= best_ok - 1e-12: break
        best_ok, cs0 = best
        cs0f = np.array(cs0, float)
    found.sort(key=lambda x: x[0])
    return found[:5]

def main():
    sin, sout, rounds, starts, seed = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
    st = json.load(open(sin))
    den = st['denominator']; num = st['numerator']
    L = Lib([dict(name=r['name'], coeffs=r['coeffs']) for r in den], [dict(name=s['name'], coeffs=s['coeffs']) for s in num], t)
    if 'mu' in st:
        mu0, w0 = np.array(st['mu']), np.array(st['w'] + [0.0]*(len(num)-len(st['w'])))
    else:
        D = sum((len(r['coeffs'])-1)*r['multiplicity'] for r in den)
        mu0 = np.array([r['multiplicity'] for r in den])/D; w0 = np.array([float(Fraction(s['q_rational'])) for s in num])/D
    W = WLP(t, L.LR, L.LM, L.dR, L.LS, L.dS, thr=0.004)
    v = W.init(mu0, w0); print('init LP', v, flush=True)
    have = {tuple(s['coeffs']) for s in num}
    dens_coeffs = [r['coeffs'] for r in den]
    rs = seed * 1000003
    for rnd in range(rounds):
        nu, price, lam = W.duals()
        wpos = np.maximum(np.array(W.w), 0) + 1e-12
        t0 = time.time()
        with Pool(int(os.environ.get('CG_PROCS', os.cpu_count())), initializer=winit,
                  initargs=(nu, dens_coeffs, [s['coeffs'] for s in num], wpos)) as pool:
            out = pool.map(irls_job, [(rs + i,) for i in range(starts)], chunksize=1)
        rs += starts
        facs = {}
        for lst in out:
            for r, cs in lst:
                for f in irr_factors(cs):
                    if tuple(f) in have: continue
                    P = fmpz_poly(f[::-1])
                    if any(P.gcd(fmpz_poly(c[::-1])).degree() > 0 for c in dens_coeffs): continue
                    facs[tuple(f)] = f
        # exact rc screen
        scored = []
        for f in facs.values():
            ls = logabs(f, t); r = nu @ ls
            if r < -1e-13: scored.append((r, f, ls))
        scored.sort(key=lambda x: x[0])
        print(f'round {rnd}: irls {time.time()-t0:.0f}s, new factors {len(facs)}, negative rc {len(scored)}', flush=True)
        added = 0
        for r, f, ls in scored[:40]:
            nu2, _, _ = W.duals(); r2 = nu2 @ ls
            if r2 >= -1e-13: continue
            W.add_num(ls, len(f)-1); v2 = W.solve()
            if v - v2 > 1e-13:
                name = f'cg{seed}_{rnd}_{added}'
                num.append(dict(name=name, coeffs=f)); have.add(tuple(f)); added += 1
                print(f'   + {name} deg {len(f)-1} rc {r2:.3e} gain {v-v2:.3e} -> {v2:.13f}', flush=True)
                v = v2
            else:
                W.pop_num(); W.solve()
        mu, w = W.weights()
        # prune: drop numerators with zero weight? keep all for now
        json.dump(dict(denominator=den, numerator=num, mu=list(mu), w=list(w), lp=v, N=N), open(sout, 'w'))
        if added == 0 and rnd > 2: pass

if __name__ == '__main__':
    main()
