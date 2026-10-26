"""Exact monochromatic-K4 density of a step graphon with equal block weights (Tao's table, C89).

Input: JSON (gz) with 'edge_probability_denominator' D and integer matrix 'red_probability_numerators' R
(symmetric, 0 <= R <= D). With Bm = D - R,
    c4 = [ sum_i tr((diag(R_i) R)^3) + sum_i tr((diag(Bm_i) Bm)^3) ] / (m^4 D^6).
The numerator N is computed modulo several primes p < 2^20.5 (float64 matrix products of residues are exact
because m p^2 < 2^53), combined by CRT, and confirmed by one extra prime. Every block is summed; no symmetry
is assumed. Written independently of the search code.

    python3 exact89.py graphon.json.gz [n_primes]
"""
import gzip
import json
import sys
import time
from fractions import Fraction

import numpy as np


def primes_below(limit, count):
    out, n = [], limit
    while len(out) < count:
        n -= 1
        if all(n % q for q in range(2, int(n ** 0.5) + 1)):
            out.append(n)
    return out


_G = {}


def _init(path):
    d = json.load(gzip.open(path))
    _G["D"] = d["edge_probability_denominator"]
    _G["R"] = np.array(d["red_probability_numerators"], dtype=np.int64)


def _rows(args):
    """partial sum over rows i in [lo, hi) of tr((D_i M)^3), for M = R and M = D - R, modulo p."""
    p, lo, hi = args
    R, D = _G["R"], _G["D"]
    tot = 0
    for M in (R % p, (D - R) % p):
        Mf = M.astype(np.float64)
        for i in range(lo, hi):
            A = np.mod(Mf * Mf[i][None, :], p)            # A_jk = M_jk M_ik  (< p)
            Q = np.mod(A @ Mf, p)                         # Q_jl = sum_k M_jk M_ik M_kl, exact: m p^2 < 2^53
            S = np.mod((Q * A).sum(axis=1), p)            # S_j = sum_l Q_jl M_jl M_il
            tot = (tot + int(np.mod(np.dot(Mf[i], S), p))) % p
    return tot


def numerator_mod(R, D, p, pool=None, path=None):
    m = R.shape[0]
    if pool is None:
        _G["R"], _G["D"] = R, D
        return _rows((p, 0, m))
    step = max(1, m // 256)
    parts = pool.map(_rows, [(p, lo, min(m, lo + step)) for lo in range(0, m, step)])
    return sum(parts) % p


def main():
    d = json.load(gzip.open(sys.argv[1]))
    D = d["edge_probability_denominator"]
    R = np.array(d["red_probability_numerators"], dtype=np.int64)
    m = R.shape[0]
    assert (R == R.T).all() and R.min() >= 0 and R.max() <= D
    bound_bits = (m ** 4 * D ** 6).bit_length() + 2
    pmax = int((2 ** 53 / m) ** 0.5)
    ps = primes_below(pmax, 64)
    need, acc = [], 1
    for p in ps:
        need.append(p)
        acc *= p
        if acc.bit_length() > bound_bits:
            break
    check = ps[len(need)]
    print(f"m={m} D={D} primes={len(need)} (each < {pmax}) + check prime {check}", flush=True)
    allp = need + [check]
    if len(sys.argv) > 2 and sys.argv[2].startswith("part="):      # part=k/K: compute only some primes
        k, K = map(int, sys.argv[2][5:].split("/"))
        import multiprocessing as mp
        workers = int(sys.argv[3]) if len(sys.argv) > 3 else 1
        pool = mp.Pool(workers, initializer=_init, initargs=(sys.argv[1],)) if workers > 1 else None
        for p in allp[k::K]:
            t0 = time.time()
            print(f"RESIDUE {p} {numerator_mod(R, D, p, pool)}  ({time.time() - t0:.0f}s)", flush=True)
        return
    if len(sys.argv) > 2 and sys.argv[2].startswith("residues="):  # combine residue lines from files
        got = {}
        for fn in sys.argv[2][9:].split(","):
            for line in open(fn):
                if line.startswith("RESIDUE"):
                    _, p, r = line.split()[:3]
                    got[int(p)] = int(r)
        res = [got[p] for p in allp]
    else:
        res = []
    t0 = time.time()
    for p in ([] if res else allp):
        r = numerator_mod(R, D, p)
        res.append(r)
        print(f"  p={p}: {r}  ({time.time() - t0:.0f}s)", flush=True)
    # CRT
    N, M = 0, 1
    for p, r in zip(need, res[:-1]):
        t = ((r - N) * pow(M, -1, p)) % p
        N, M = N + M * t, M * p
    assert N % check == res[-1], "check prime disagrees"
    c4 = Fraction(N, m ** 4 * D ** 6)
    print(f"EXACT c4 = {c4.numerator}/{c4.denominator}")
    print(f"        = {float(c4):.16f}")


if __name__ == "__main__":
    main()
