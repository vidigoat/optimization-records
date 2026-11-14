"""Warm-start LP (highspy) for column generation on the C82 Doche LP."""
import numpy as np, highspy
from core import logabs, logmahlerT, roots

class WLP:
    def __init__(self, t, LR, LM, dR, LS, dS, thr=0.004, budget=1-1e-7):
        self.t, self.N = t, len(t)
        self.LR, self.LM, self.dR = LR, LM, dR
        self.LS, self.dS = [x for x in LS], list(dS)
        self.thr, self.budget = thr, budget
        self.h = highspy.Highs(); self.h.setOptionValue('output_flag', False)
        self.h.setOptionValue('primal_feasibility_tolerance', 1e-10)
        self.h.setOptionValue('dual_feasibility_tolerance', 1e-10)
        self.rows = []   # grid indices of epigraph rows
    def build(self, act):
        h = self.h; h.clearModel(); inf = highspy.kHighsInf
        J, I = len(self.LM), len(self.LS)
        # rows: 0 eq (dR.mu=1), 1 degree budget, then epigraph rows
        h.addRows(2, np.array([1.0, -inf]), np.array([1.0, self.budget]), 0, np.array([0,0], dtype=np.int32), np.array([], dtype=np.int32), np.array([]))
        self.rows = list(act)
        n = len(act)
        h.addRows(n, np.full(n, -inf), np.zeros(n), 0, np.zeros(n, dtype=np.int32), np.array([], dtype=np.int32), np.array([]))
        # mu columns
        for j in range(J):
            idx = np.concatenate([[0], 2+np.arange(n)]).astype(np.int32)
            val = np.concatenate([[self.dR[j]], -self.LR[j][act]])
            h.addCol(self.LM[j], 0, inf, len(idx), idx, val)
        for i in range(I):
            idx = np.concatenate([[1], 2+np.arange(n)]).astype(np.int32)
            val = np.concatenate([[self.dS[i]], self.LS[i][act]])
            h.addCol(0.0, 0, inf, len(idx), idx, val)
        for k in range(n):
            h.addCol(1.0/self.N, 0, inf, 1, np.array([2+k], dtype=np.int32), np.array([-1.0]))
        self.J, self.I = J, I
    def col_order(self):
        # column indices: mu 0..J-1, w J..J+I-1 (appended numerators are placed after y's -> track)
        pass
    def solve(self, maxit=20):
        h = self.h
        for it in range(maxit):
            h.run()
            st = h.getModelStatus()
            assert st == highspy.HighsModelStatus.kOptimal, st
            mu, w = self.weights()
            g = self.g(mu, w)
            have = np.zeros(self.N, bool); have[self.rows] = True
            new = np.flatnonzero((g > -self.thr/4) & ~have)
            if len(new) == 0: break
            self.add_rows(np.flatnonzero((g > -self.thr) & ~have))
        self.mu, self.w = mu, w
        return h.getInfo().objective_function_value
    def add_rows(self, idx):
        h = self.h; inf = highspy.kHighsInf
        n0 = h.getNumRow()
        cols_mu, cols_w = self.colmap()
        for k, t in enumerate(idx):
            ind = list(cols_mu) + list(cols_w); val = [-x[t] for x in self.LR] + [x[t] for x in self.LS]
            h.addRow(-inf, 0.0, len(ind), np.array(ind, dtype=np.int32), np.array(val))
        for k in range(len(idx)):
            h.addCol(1.0/self.N, 0, inf, 1, np.array([n0+k], dtype=np.int32), np.array([-1.0]))
        self.rows += list(idx)
    def colmap(self):
        # mu columns 0..J-1; w columns: first I0 at J..J+I0-1, the rest recorded in self.extra_cols
        return self._mu_cols, self._w_cols
    def g(self, mu, w):
        return np.dot(w, np.array(self.LS)) - np.dot(mu, self.LR) if len(self.LS) else -np.dot(mu, self.LR)
    def weights(self):
        x = np.array(self.h.getSolution().col_value)
        return x[self._mu_cols], x[self._w_cols]
    def init(self, mu0=None, w0=None):
        if mu0 is None:
            act = np.arange(self.N)
        else:
            act = np.flatnonzero(self.g(mu0, w0) > -self.thr)
        self.build(act)
        J, I, n = self.J, self.I, len(act)
        self._mu_cols = list(range(J)); self._w_cols = list(range(J, J+I))
        return self.solve()
    def add_num(self, ls, deg):
        h = self.h; inf = highspy.kHighsInf
        rows = np.array(self.rows)
        idx = np.concatenate([[1], 2+np.arange(len(rows))]).astype(np.int32)
        val = np.concatenate([[deg], ls[rows]])
        h.addCol(0.0, 0, inf, len(idx), idx, val)
        self._w_cols.append(h.getNumCol()-1); self.LS.append(ls); self.dS.append(deg)
    def pop_num(self):
        c = self._w_cols.pop(); self.LS.pop(); self.dS.pop()
        self.h.deleteCols(1, np.array([c], dtype=np.int32))
        # columns after c shift down by one
        self._w_cols = [x-1 if x > c else x for x in self._w_cols]
        self._mu_cols = [x-1 if x > c else x for x in self._mu_cols]
    def duals(self):
        sol = self.h.getSolution(); rd = np.array(sol.row_dual)
        nu = np.zeros(self.N); nu[self.rows] = -rd[2:]
        return nu, -rd[1], rd[0]

def _add_den(self, lr, lm, deg):
    h = self.h; inf = highspy.kHighsInf
    rows = np.array(self.rows)
    idx = np.concatenate([[0], 2+np.arange(len(rows))]).astype(np.int32)
    val = np.concatenate([[deg], -lr[rows]])
    h.addCol(lm, 0, inf, len(idx), idx, val)
    self._mu_cols.append(h.getNumCol()-1); self.LR = np.vstack([self.LR, lr[None, :]])
    self.LM = np.append(self.LM, lm); self.dR = np.append(self.dR, deg)

def _pop_den(self):
    c = self._mu_cols.pop(); self.LR = self.LR[:-1]; self.LM = self.LM[:-1]; self.dR = self.dR[:-1]
    self.h.deleteCols(1, np.array([c], dtype=np.int32))
    self._w_cols = [x-1 if x > c else x for x in self._w_cols]
    self._mu_cols = [x-1 if x > c else x for x in self._mu_cols]

WLP.add_den = _add_den
WLP.pop_den = _pop_den
