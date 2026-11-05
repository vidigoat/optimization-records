"""Dense primal-dual (Mehrotra) IPM for the C86 LP, same model and API as fulllp.FullLP.

    min c.w   s.t.  G w >= h  (potential cuts, integral rows, degree row as -sum b deg >= -(1-1e-8)),
                    sum(tri weights) = 1,  w >= 0.
Normal equations K = G^T diag(y/s) G + diag(z/x) are formed and factored densely (Accelerate BLAS/LAPACK).
"""
import os, time
import numpy as np
import scipy.linalg as sla
from pl import H


def tri_pot(l, c, r, z):
    l, c, r = map(np.asarray, (l, c, r)); h = 2 / (r - l)
    s = [h / (c - l), -h / (c - l) - h / (r - c), h / (r - c)]
    z = np.asarray(z)
    out = np.zeros((len(z), len(l)))
    if not np.iscomplexobj(z):
        # real z: Re H(w) = w^2 (log|w|/2 - 3/4) (the i*pi of log w for w < 0 drops out of the real part)
        for sj, kj in zip(s, (l, c, r)):
            w = z[:, None] - kj[None, :]
            aw = np.abs(w); lw = np.log(np.where(aw > 0, aw, 1.0))
            out += sj[None, :] * (w * w * (0.5 * lw - 0.75))
        return out
    for sj, kj in zip(s, (l, c, r)):
        out += (sj[None, :] * H(z[:, None] - kj[None, :])).real
    return out


def logabs(Q, x):
    return np.log(np.abs(np.polynomial.polynomial.polyval(np.asarray(x, float), np.asarray(Q, float))))


def tri_int(Q, l, c, r):
    Q = np.asarray(Q, float)
    roots = np.polynomial.polynomial.polyroots(Q)
    return np.log(abs(Q[-1])) + tri_pot(l, c, r, roots).sum(axis=0)


def ipm(G, h, A, b, c, tol=1e-10, maxit=200, verbose=False):
    """min c x s.t. G x >= h, A x = b, x >= 0. Returns x, y (>=0, for G rows), lam (for A rows), info."""
    m, n = G.shape; me = A.shape[0]
    x = np.ones(n); s = np.maximum(G @ x - h, 1.0); y = np.ones(m); z = np.ones(n); lam = np.zeros(me)
    nb = 1 + max(np.abs(h).max(), np.abs(b).max()); nc = 1 + np.abs(c).max()
    for it in range(maxit):
        rd = c - G.T @ y - A.T @ lam - z
        rp = G @ x - s - h
        re = A @ x - b
        mu = (s @ y + x @ z) / (m + n)
        pobj = c @ x; dobj = h @ y + b @ lam
        gap = abs(pobj - dobj) / (1 + abs(pobj))
        pres = max(np.abs(rp).max() / nb, np.abs(re).max() / nb); dres = np.abs(rd).max() / nc; res = max(pres, dres)
        if verbose: print(f'  it {it:3d} pobj {pobj:.12f} dobj {dobj:.12f} gap {gap:.2e} pres {pres:.2e} dres {dres:.2e} mu {mu:.2e}', flush=True)
        if gap < tol and res < tol: break
        if mu < 1e-14 * (1 + abs(pobj)) and pres < 1e-9: break
        W = y / s; D = z / x
        K = (G.T * W) @ G
        K[np.diag_indices(n)] += D + 1e-14 * (1 + np.abs(K.diagonal()))
        try:
            cf = sla.cho_factor(K, lower=True, check_finite=False)
        except np.linalg.LinAlgError:
            K[np.diag_indices(n)] += 1e-10 * (1 + np.abs(K.diagonal())); cf = sla.cho_factor(K, lower=True, check_finite=False)
        KA = sla.cho_solve(cf, A.T, check_finite=False); S = A @ KA

        def direction(rys, rzx):
            q = -rd - G.T @ ((rys + y * rp) / s) - rzx / x
            Kq = sla.cho_solve(cf, q, check_finite=False)
            dl = np.linalg.solve(S, -re - A @ Kq)
            dx = Kq + KA @ dl
            ds = G @ dx + rp
            dy = -(rys + y * ds) / s
            dz = -(rzx + z * dx) / x
            return dx, ds, dy, dz, dl

        def step(v, dv):
            neg = dv < 0
            return min(1.0, (-v[neg] / dv[neg]).min()) if neg.any() else 1.0
        dx, ds, dy, dz, dl = direction(s * y, x * z)
        ap = min(step(x, dx), step(s, ds)); ad = min(step(y, dy), step(z, dz))
        mua = ((s + ap * ds) @ (y + ad * dy) + (x + ap * dx) @ (z + ad * dz)) / (m + n)
        sig = (mua / mu) ** 3
        dx2, ds2, dy2, dz2, dl2 = direction(s * y + ds * dy - sig * mu, x * z + dx * dz - sig * mu)
        ap = 0.995 * min(step(x, dx2), step(s, ds2)); ad = 0.995 * min(step(y, dy2), step(z, dz2))
        ap = min(ap, 1.0); ad = min(ad, 1.0)
        x += ap * dx2; s += ap * ds2; y += ad * dy2; z += ad * dz2; lam += ad * dl2
    return x, y, lam, dict(it=it, gap=gap, res=res, pobj=pobj, dobj=dobj)


class DenseLP:
    def __init__(self, margin=2e-6, imargin=1e-7, verbose=False):
        self.margin, self.imargin, self.verbose = margin, imargin, verbose
        self.tris = np.zeros((0, 3)); self.polys = []; self.x = np.zeros(0)
        self.Ut = np.zeros((0, 0))       # potentials of tris at points (np, nt)
        self.Lq = np.zeros((0, 0))       # -log|Q| at points (np, nq)
        self.It = np.zeros((0, 0))       # integrals (nq, nt)

    def add_tris(self, T):
        T = np.atleast_2d(np.asarray(T, float))
        if len(T) == 0: return
        l, c, r = T.T
        self.Ut = np.hstack([self.Ut.reshape(len(self.x), len(self.tris)), tri_pot(l, c, r, self.x)])
        if self.polys:
            self.It = np.hstack([self.It, np.array([tri_int(Q, l, c, r) for Q in self.polys])])
        else:
            self.It = np.zeros((0, len(self.tris) + len(T)))
        self.tris = np.vstack([self.tris, T])

    def add_polys(self, Qs):
        for Q in Qs:
            Q = [int(t) for t in Q]
            T = self.tris
            self.It = np.vstack([self.It.reshape(len(self.polys), len(T)), tri_int(Q, T[:, 0], T[:, 1], T[:, 2])[None, :]])
            self.Lq = np.hstack([self.Lq.reshape(len(self.x), len(self.polys)), -logabs(Q, self.x)[:, None]])
            self.polys.append(Q)

    def add_points(self, x):
        x = np.asarray(x, float)
        if len(x) == 0: return
        T = self.tris
        Ut = tri_pot(T[:, 0], T[:, 1], T[:, 2], x) if len(T) else np.zeros((len(x), 0))
        Lq = np.column_stack([-logabs(Q, x) for Q in self.polys]) if self.polys else np.zeros((len(x), 0))
        self.Ut = np.vstack([self.Ut.reshape(len(self.x), len(T)), Ut]); self.Lq = np.vstack([self.Lq.reshape(len(self.x), len(self.polys)), Lq])
        self.x = np.concatenate([self.x, x])

    def solve(self, tol=1e-10):
        nt, nq, npt = len(self.tris), len(self.polys), len(self.x)
        deg = np.array([len(Q) - 1 for Q in self.polys], float)
        G = np.zeros((npt + nq + 1, nt + nq))
        G[:npt, :nt] = self.Ut; G[:npt, nt:] = self.Lq
        G[npt:npt + nq, :nt] = self.It
        G[-1, nt:] = -deg
        h = np.concatenate([self.margin * np.ones(npt), self.imargin * np.ones(nq), [-(1 - 1e-8)]])
        A = np.zeros((1, nt + nq)); A[0, :nt] = 1
        c = np.concatenate([self.tris.mean(axis=1), np.zeros(nq)])
        t = time.time()
        w, y, lam, info = ipm(G, h, A, np.array([1.0]), c, tol=tol, verbose=self.verbose)
        self.solve_time = time.time() - t; self.info = info
        self.w_tri = w[:nt]; self.b = w[nt:]
        self.y_pot = y[:npt]; self.y_int = y[npt:npt + nq]; self.y_deg = y[-1]; self.y_mass = lam[0]
        self.cost = c
        return info['pobj']

    def f(self, x):
        x = np.asarray(x, float); out = np.zeros(len(x))
        # slope jumps of the piecewise-linear density on the union of knots: U(x) = sum_k J_k Re H(x - e_k)
        T = self.tris; w = self.w_tri; h = 2 / (T[:, 2] - T[:, 0])
        js = np.concatenate([w * h / (T[:, 1] - T[:, 0]), -w * h / (T[:, 1] - T[:, 0]) - w * h / (T[:, 2] - T[:, 1]), w * h / (T[:, 2] - T[:, 1])])
        ks = np.concatenate([T[:, 0], T[:, 1], T[:, 2]])
        kn, inv = np.unique(ks, return_inverse=True); J = np.bincount(inv, weights=js)
        def chunk(i):
            wv = x[i:i + 2000, None] - kn[None, :]; aw = np.abs(wv)
            out[i:i + 2000] = (wv * wv * (0.5 * np.log(np.where(aw > 0, aw, 1.0)) - 0.75)) @ J
        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor(int(os.environ.get('FTHREADS', 12))) as ex: list(ex.map(chunk, range(0, len(x), 2000)))
        for Q, bj in zip(self.polys, self.b):
            out -= bj * logabs(Q, x)
        return out

    def tri_reduced_costs(self, T):
        T = np.atleast_2d(np.asarray(T, float)); l, c, r = T[:, 0], T[:, 1], T[:, 2]
        rcost = (l + c + r) / 3 - self.y_mass - tri_pot(l, c, r, self.x).T @ self.y_pot
        for Q, y in zip(self.polys, self.y_int):
            rcost -= y * tri_int(Q, l, c, r)
        return rcost
