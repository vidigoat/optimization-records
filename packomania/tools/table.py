"""Packomania's published records, read from the files downloaded into runs/packomania."""
import os

ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "runs", "packomania")


def radii(kind):
    out = {}
    with open(os.path.join(ROOT, f"{kind}_radius.txt"), errors="replace") as f:
        for line in f:
            p = line.split()
            if len(p) >= 2 and p[0].isdigit():
                out[int(p[0])] = p[1]
    return out


def authors(kind):
    out = {}
    with open(os.path.join(ROOT, f"{kind}_author.txt"), errors="replace") as f:
        for line in f:
            p = line.split(None, 1)
            if len(p) == 2 and p[0].isdigit():
                out[int(p[0])] = p[1].strip()
    return out


def coords(kind, n):
    """The published packing for n, as tuples of decimal strings (the trailing radius column dropped)."""
    from .verify import DIM
    d, pts = DIM[kind], []
    with open(os.path.join(ROOT, f"{kind}_coords", f"{kind}{n}.txt")) as f:
        for line in f:
            if line.lstrip().startswith("#"):
                continue
            p = line.split()
            if len(p) >= d + 1:
                pts.append(tuple(p[1:1 + d]))
    return pts


def dominated(kind):
    """Records beaten outright by a larger packing: {n: (donor m, donor radius)} where m > n and r(m) > r(n).

    Deleting m - n objects from a valid packing of m leaves a valid packing of n with the same radius.
    """
    rad = {n: float(r) for n, r in radii(kind).items()}
    have = {n for n in rad if os.path.exists(os.path.join(ROOT, f"{kind}_coords", f"{kind}{n}.txt"))}
    out, best, best_n = {}, 0.0, None
    for n in sorted(rad, reverse=True):
        if rad[n] < best:
            out[n] = (best_n, best)
        if rad[n] > best and n in have:          # a donor is only useful if its coordinates are published
            best, best_n = rad[n], n
    return out
