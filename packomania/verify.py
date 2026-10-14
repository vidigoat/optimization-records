"""The trusted checker for packings, in exact arithmetic.

Containers follow Packomania's conventions: a square of side 1 centred at the origin ("csq"), a
circle of radius 1 ("cci"), a cube of edge 1 ("scu"), a sphere of radius 1 ("ssp"). A packing of
N equal objects of radius r is valid when every object lies inside the container and no two
overlap. Coordinates are read as exact decimals (Fractions), so rounding can never let an
overlapping packing through: touching is allowed, overlapping by any amount is not.
"""
import math
from fractions import Fraction

DIM = {"csq": 2, "cci": 2, "scu": 3, "ssp": 3}


def exact(x):
    return x if isinstance(x, Fraction) else Fraction(str(x))


def inside(kind, p, r):
    if kind in ("csq", "scu"):
        half = Fraction(1, 2)
        return all(abs(c) <= half - r for c in p)
    return r <= 1 and sum(c * c for c in p) <= (1 - r) ** 2        # cci, ssp


def verify(kind, points, radius):
    """Returns (ok, reason). `points` are tuples of numbers or decimal strings; `radius` likewise."""
    d = DIM[kind]
    r = exact(radius)
    if r <= 0:
        return False, "radius must be positive"
    pts = [tuple(exact(c) for c in p) for p in points]
    if any(len(p) != d for p in pts):
        return False, f"every point needs {d} coordinates"
    for i, p in enumerate(pts):
        if not inside(kind, p, r):
            return False, f"object {i + 1} sticks out of the container"
    # Only nearby pairs can overlap: bucket by a grid of cell size just over 2r (floats only choose pairs).
    cell = 2 * float(r) * 1.000001
    grid = {}
    for i, p in enumerate(pts):
        grid.setdefault(tuple(math.floor(float(c) / cell) for c in p), []).append(i)
    four_r2 = 4 * r * r
    offsets = [()]
    for _ in range(d):
        offsets = [o + (k,) for o in offsets for k in (-1, 0, 1)]
    for key, members in grid.items():
        for off in offsets:
            other = grid.get(tuple(a + b for a, b in zip(key, off)))
            if not other:
                continue
            for i in members:
                for j in other:
                    if i < j and sum((a - b) ** 2 for a, b in zip(pts[i], pts[j])) < four_r2:
                        return False, f"objects {i + 1} and {j + 1} overlap"
    return True, f"{len(pts)} objects of radius {float(r):.12f} fit"


def largest_radius(kind, points):
    """The largest radius the given centres allow, in floats (for choosing a radius to then verify exactly)."""
    pts = [tuple(float(c) for c in p) for p in points]
    if kind in ("csq", "scu"):
        wall = min(0.5 - abs(c) for p in pts for c in p)
    else:
        wall = min(1 - math.sqrt(sum(c * c for c in p)) for p in pts)
    cell = 2 * wall + 1e-9
    grid = {}
    for i, p in enumerate(pts):
        grid.setdefault(tuple(math.floor(c / cell) for c in p), []).append(i)
    d = len(pts[0])
    offsets = [()]
    for _ in range(d):
        offsets = [o + (k,) for o in offsets for k in (-1, 0, 1)]
    half = wall
    for key, members in grid.items():
        for off in offsets:
            for j in grid.get(tuple(a + b for a, b in zip(key, off)), ()):
                for i in members:
                    if i < j:
                        half = min(half, math.dist(pts[i], pts[j]) / 2)
    return min(wall, half)
