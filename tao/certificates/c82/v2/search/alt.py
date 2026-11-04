"""Alternate numerator CG (cg.irls_job) and denominator CG (dsearch.job2) on one warm LP.
usage: alt.py state_in state_out rounds num_starts den_starts seed"""
import numpy as np, json, sys, time, os
from multiprocessing import Pool
from flint import fmpz_poly
import cg, dsearch
from core import *
from wlp import WLP
t = cg.t; N = cg.N
def init_both(nu, lam, dens, lib, wpos, spool, nump, near):
    cg.winit(nu, dens, lib, wpos)
    dsearch.winit2(nu, lam, spool, nump, dens, near)
def numjob(s):
    return cg.irls_job((s,))
def main():
    sin, sout, rounds, ns, nd, seed = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
    st = json.load(open(sin)); den, num = st['denominator'], st['numerator']
    L = Lib([dict(name=r['name'], coeffs=r['coeffs']) for r in den], [dict(name=s['name'], coeffs=s['coeffs']) for s in num], t)
    W = WLP(t, L.LR, L.LM, L.dR, L.LS, L.dS, thr=0.004)
    w0 = np.array(st['w'] + [0.0]*(len(num)-len(st['w']))); mu0 = np.array(st['mu'] + [0.0]*(len(den)-len(st['mu'])))
    v = W.init(mu0, w0); print('init LP', v, flush=True)
    have = {tuple(r['coeffs']) for r in den} | {tuple(s['coeffs']) for s in num}
    rs = seed*1000003; procs = int(os.environ.get('CG_PROCS', os.cpu_count()))
    for rnd in range(rounds):
        for phase in ('num', 'den'):
            nu, price, lam = W.duals(); mu, w = W.weights()
            dens = [r['coeffs'] for r in den]
            wpos = np.maximum(np.array(w), 0) + 1e-12
            act = [r['coeffs'] for r, x in zip(den, mu) if x > 1e-9]
            spool = act + [r['coeffs'] for r in den if len(r['coeffs']) - 1 <= 28]
            near = [r['coeffs'] for r, lr_, lm_ in zip(den, W.LR, W.LM) if lm_ - nu @ lr_ - lam*(len(r['coeffs'])-1) < float(os.environ.get('DS_NEARTHR', 0.003))]
            n = ns if phase == 'num' else nd
            if n == 0: continue
            t0 = time.time()
            with Pool(procs, initializer=init_both, initargs=(nu, lam, dens, [s['coeffs'] for s in num], wpos, spool, [s['coeffs'] for s in num], near)) as pool:
                out = pool.map(numjob if phase == 'num' else dsearch.job2, [rs+i for i in range(n)], chunksize=1)
            rs += n
            facs = {}
            for lst in out:
                for item in lst:
                    if item[0] == 'best': continue
                    r, cs = item
                    for f, e in fmpz_poly(cs[::-1]).factor()[1]:
                        if f.degree() < (1 if phase == 'num' else 2): continue
                        c = [int(x) for x in f.coeffs()][::-1]
                        if c[0] < 0: c = [-x for x in c]
                        if tuple(c) in have: continue
                        if phase == 'den' and abs(c[0]) != 1: continue
                        if phase == 'num' and any(fmpz_poly(c[::-1]).gcd(fmpz_poly(q[::-1])).degree() > 0 for q in dens): continue
                        facs[tuple(c)] = c
            scored = []
            for c in facs.values():
                lr = logabs(c, t)
                if phase == 'num': r = nu @ lr + price*(len(c)-1); lm = None
                else: lm = logmahlerT(c); r = lm - nu @ lr - lam*(len(c)-1)
                if r < -1e-12: scored.append((r, c, lr, lm))
            scored.sort(key=lambda x: x[0])
            bests = [o[-1][1] for o in out if phase == 'den' and o and o[-1][0] == 'best' and o[-1][1] is not None]
            print(f'round {rnd} {phase}: {time.time()-t0:.0f}s starts {n} cands {len(facs)} neg {len(scored)}' + (f' near {len(near)} best-rc q ' + str(np.round(np.quantile(bests, [0, .1, .5]), 5)) if bests else ''), flush=True)
            added = 0
            batch = []
            nu2, price2, lam2 = W.duals()
            for r, c, lr, lm in scored[:int(os.environ.get('ALT_BATCH', 60))]:
                if phase == 'num': W.add_num(lr, len(c)-1)
                else: W.add_den(lr, lm, len(c)-1)
                batch.append((r, c, lm))
            if batch:
                v2 = W.solve(); mu_, w_ = W.weights()
                wts = (w_ if phase == 'num' else mu_)[-len(batch):]
                for (r, c, lm), x in zip(batch, wts):
                    name = ('n' if phase == 'num' else 'D') + f'{seed}_{rnd}_{added}'; added += 1
                    (num if phase == 'num' else den).append(dict(name=name, coeffs=c)); have.add(tuple(c))
                    if x > 1e-12: print(f'   + {name} deg {len(c)-1} rc {r:.3e} weight {x:.3e}' + (f' logM/deg {lm/(len(c)-1):.9f}' if lm else ''), flush=True)
                print(f'   batch {len(batch)} gain {v-v2:.3e} -> {v2:.13f}', flush=True)
                v = v2
            mu, w = W.weights()
            json.dump(dict(denominator=den, numerator=num, mu=list(mu), w=list(w), lp=v, N=N), open(sout, 'w'))
if __name__ == '__main__':
    main()
