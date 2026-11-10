"""C86 search with the dense IPM: Pinho's node grid as triangle columns, cutting planes, triangle and auxiliary pricing.
Checkpoints after every solve.   python run_dense.py out_prefix [extra_polys.json ...]"""
import sys, os, json, time, pickle
import numpy as np
from fractions import Fraction
from multiprocessing import Pool
sys.path.insert(0, '.')
from denselp import DenseLP, logabs
from auxsearch import search, rc

G = {}
def winit(x, y, z, lib, libw): G.update(x=x, y=y, z=z, lib=lib, libw=libw)
def job(seed): return search(G['x'], G['y'], G['z'], 3, np.random.default_rng(seed), G['lib'], G['libw'], dmin=8, dmax=int(os.environ.get('DMAX', 30)))


def refine(L, xx, fx, thr=float(os.environ.get('RTHR', 1e-5))):
    """Locate every grid local minimum with f < thr to ~1e-12 by repeated sub-gridding of its bracket."""
    i = np.nonzero((fx[1:-1] <= fx[:-2]) & (fx[1:-1] <= fx[2:]) & (fx[1:-1] < thr))[0] + 1
    if len(i) == 0: return np.zeros(0), np.zeros(0)
    a, b = xx[i - 1], xx[i + 1]
    for _ in range(5):
        sub = a[:, None] + (b - a)[:, None] * np.linspace(0, 1, 21)[None, :]
        fs = L.f(sub.ravel()).reshape(sub.shape)
        j = fs.argmin(axis=1); r = np.arange(len(a))
        a, b = sub[r, np.maximum(j - 1, 0)], sub[r, np.minimum(j + 1, 20)]
    return sub[r, j], fs[r, j]


def save(L, out, v, tag):
    pickle.dump(dict(tris=L.tris, w=L.w_tri, polys=L.polys, b=L.b, x=L.x, v=v), open(f'{out}_state.pkl', 'wb'))
    print(f'   [{tag}] LP {v:.10f} ({L.solve_time:.0f}s, {L.info["it"]} it)  {len(L.tris)} tris {len(L.polys)} polys {len(L.x)} pts', flush=True)


def main():
    out = sys.argv[1]
    LO, HI = 0.04, 5.30
    if os.environ.get('RESUME'):
        st = pickle.load(open(os.environ['RESUME'], 'rb'))
        L = DenseLP(); L.add_points(st['x']); 
        k = st['w'] > 1e-12 * st['w'].max(); L.add_tris(st['tris'][k] if os.environ.get('PRUNE') else st['tris'])
        L.add_polys(st['polys'])
    else:
        c = json.load(open('../pinho/certificate.json'))
        e = np.array([float(Fraction(v)) for v in c['edges']])
        k = np.arange(1, len(e) - 1); sel = (e[k] > LO) & (e[k] < HI)
        T = np.column_stack([e[k - 1], e[k], e[k + 1]])[sel]
        pts = np.unique(np.concatenate([e[(e > LO - 0.01) & (e < HI + 0.01)], ((e[1:] + e[:-1]) / 2)[(e[1:] > LO - 0.01) & (e[:-1] < HI + 0.01)]]))
        L = DenseLP(); L.add_points(pts); L.add_tris(T); L.add_polys([list(a) for a in c['A']])
    for fn in sys.argv[2:]:
        have = {tuple(q) for q in L.polys}
        L.add_polys([p for p in json.load(open(fn)) if tuple(p) not in have])
    v = L.solve(); save(L, out, v, 'init')

    LIB = []
    libdir = os.environ.get('LIBDIR', '/Users/vidigoat/Library/Application Support/Claude/scratch-workspaces/74f4ab68-8952-4986-a614-17afc140842a/e7c77990-ad0d-466c-8fce-b86019184e88/scratch-2026-10-13-38ad35/optimization-records/tao/certificates/c86')
    for fn in ['deg1to7.txt', 'deg8.txt']:
        for line in open(libdir + '/' + fn):
            if line.strip(): LIB.append(json.loads(line))
    for fn in os.environ.get('XLIB', '').split(','):
        if fn:
            for line in open(fn):
                if line.strip().startswith('['): LIB.append(json.loads(line))
    seed = int(os.environ.get('SEED', 1000))
    xx = np.linspace(float(os.environ.get('XLO', 0.015)), float(os.environ.get('XHI', 5.65)), int(os.environ.get('NX', 1000001)))
    for rnd in range(int(os.environ.get('ROUNDS', 200))):
        t0 = time.time()
        for _ in range(12):
            fx = L.f(xx)
            t1 = time.time(); xm, fm = refine(L, xx, fx); print(f'   refine: {len(xm)} minima in {time.time() - t1:.0f}s', flush=True)
            idx = fm < L.margin / 2
            print(f'   cut: min f {fx.min():.2e} (refined {fm.min() if len(fm) else 0:.2e}), {idx.sum()} refined local minima below margin/2', flush=True)
            if idx.sum() == 0: break
            L.add_points(xm[idx]); v = L.solve(); save(L, out, v, 'cut')
        fmin = L.f(xx).min()
        knots = np.unique(L.tris[L.w_tri > 1e-12 * L.w_tri.max()].ravel())
        knots = knots[(knots > LO) & (knots < HI)]
        cand = []
        for a, b in zip(knots[:-1], knots[1:]):
            m = (a + b) / 2
            cand += [(a, m, b), (a, (a + m) / 2, m), (m, (m + b) / 2, b)]
        cand = np.array(cand)
        rcs = np.concatenate([L.tri_reduced_costs(cand[i:i + 2000]) for i in range(0, len(cand), 2000)])
        neg = np.argsort(rcs)[:int(os.environ.get('NTRI', 300))]; neg = neg[rcs[neg] < -1e-9]
        y = np.maximum(L.y_pot, 0); kk = y > 1e-12 * y.max()
        have = {tuple(q) for q in L.polys}
        lr = sorted([(rc(f, L.x, y, L.y_deg), f) for f in LIB if tuple(f) not in have])
        newp = [f for r_, f in lr[:20] if r_ < -1e-9]
        with Pool(int(os.environ.get('NPROC', 8)), initializer=winit, initargs=(L.x[kk], y[kk], L.y_deg, [list(q) for q in L.polys], list(np.maximum(L.b, 0)))) as pool:
            nj = int(os.environ.get('NJOBS', 16)); res = pool.map(job, range(seed, seed + nj)); seed += nj
        found = {}
        for lst in res:
            for r_, f in lst:
                if tuple(f) not in have: found[tuple(f)] = r_
        newp += [list(f) for f, r_ in sorted(found.items(), key=lambda t: t[1])[:int(os.environ.get('NPOLY', 15))]]
        if len(neg): L.add_tris(cand[neg])
        seen = {tuple(q) for q in L.polys}; uniq = []
        for q in newp:
            q = [int(t) for t in q]
            if q[-1] < 0: q = [-t for t in q]
            if tuple(q) not in seen: seen.add(tuple(q)); uniq.append(q)
        newp = uniq
        if newp: L.add_polys(newp)
        if len(neg) or newp:
            v = L.solve(); save(L, out, v, 'price')
        json.dump(L.polys, open(f'{out}_polys.json', 'w'))
        print(f'round {rnd}: LP {v:.10f}  min f {fmin:.2e}  +{len(neg)} triangles (best rc {rcs.min():.2e})  +{len(newp)} polys '
              f'(best rc {min([r for r, _ in lr[:1]] + list(found.values()) + [0]):.2e})  {time.time() - t0:.0f}s', flush=True)
        if not len(neg) and not newp: break


if __name__ == '__main__':
    main()
