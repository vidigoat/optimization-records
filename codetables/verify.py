"""Trusted checker for linear codes over GF(q), q in {2, 3, 4, 5, 7, 8}.

A code is given by a k x n generator matrix with entries 0..q-1. For q = 4 and 8 an entry is a bit
vector in the polynomial basis (GF(4): alpha^2 = alpha + 1; GF(8): alpha^3 = alpha + 1). The checker
confirms the rank is k, then computes the minimum distance exactly by enumerating every codeword
(bin/mindist, built from codes/mindist.c). Nothing is taken from the program that proposed the code.
"""
import os
import subprocess

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(ROOT, "bin", "mindist")
FIELDS = (2, 3, 4, 5, 7, 8)


def mul(q, a, b):
    if q in (2, 3, 5, 7):
        return a * b % q
    m, poly, r = (2, 0b111, 0) if q == 4 else (3, 0b1011, 0)
    for i in range(m):
        if b >> i & 1:
            r ^= a << i
    for i in range(2 * m - 2, m - 1, -1):
        if r >> i & 1:
            r ^= poly << (i - m)
    return r


def add(q, a, b):
    return a ^ b if q in (2, 4, 8) else (a + b) % q


def inv(q, a):
    return next(x for x in range(1, q) if mul(q, a, x) == 1)


def neg(q, a):
    return a if q in (2, 4, 8) else (-a) % q


def rank(q, G):
    M = [row[:] for row in G]
    r = 0
    cols = len(M[0]) if M else 0
    for c in range(cols):
        p = next((i for i in range(r, len(M)) if M[i][c]), None)
        if p is None:
            continue
        M[r], M[p] = M[p], M[r]
        f = inv(q, M[r][c])
        M[r] = [mul(q, f, x) for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c]:
                g = neg(q, M[i][c])
                M[i] = [add(q, x, mul(q, g, y)) for x, y in zip(M[i], M[r])]
        r += 1
    return r


def min_distance(q, G, target=0, timeout=3600):
    k, n = len(G), len(G[0])
    data = f"{q} {n} {k} {target}\n" + " ".join(str(x) for row in G for x in row) + "\n"
    out = subprocess.run([BIN], input=data, capture_output=True, text=True, timeout=timeout)
    if out.returncode != 0:
        raise ValueError(f"checker failed ({out.returncode})")
    d, count = out.stdout.split()
    return int(d), int(count)


def check(q, G, claimed_d=None):
    """Returns (ok, d, reason). ok means: valid entries, full rank, and (if given) distance >= claimed_d."""
    if q not in FIELDS:
        return False, 0, f"q must be one of {FIELDS}"
    if not G or any(len(r) != len(G[0]) for r in G):
        return False, 0, "generator matrix must be rectangular"
    if any(not isinstance(x, int) or not 0 <= x < q for r in G for x in r):
        return False, 0, "entries must be integers 0..q-1"
    k, n = len(G), len(G[0])
    if rank(q, G) != k:
        return False, 0, "rows are not linearly independent"
    if q ** k > 4e9:
        return False, 0, "too many codewords to check exhaustively"
    d, _ = min_distance(q, G, target=claimed_d or 0)
    if claimed_d is not None and d < claimed_d:
        return False, d, f"has a codeword of weight {d} < {claimed_d}"
    return True, d, f"[{n},{k},{d}] over GF({q})"
