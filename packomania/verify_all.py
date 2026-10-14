"""Check every packing in a submission zip against the published records, in exact arithmetic.

Each <kind><n>.pck file holds the radius (line 1), the author (line 2) and one line of coordinates per
object. SUMMARY.tsv in the zip gives the published radius for each entry at the time of submission.

    python3 packomania/verify_all.py packomania/packings.zip
"""
import os
import sys
import zipfile
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from verify import DIM, verify  # noqa: E402


def main(path):
    z = zipfile.ZipFile(path)
    published = {}
    for line in z.read("SUMMARY.tsv").decode().splitlines()[1:]:
        kind, n, old, new, gain = line.split("\t")
        published[(kind, int(n))] = old
    good = bad = 0
    for name in sorted(z.namelist()):
        if not name.endswith(".pck"):
            continue
        stem = name[:-4]
        kind = next(k for k in DIM if stem.startswith(k) and stem[len(k):].isdigit())
        n = int(stem[len(kind):])
        lines = z.read(name).decode().splitlines()
        radius, points = lines[0].strip(), [tuple(l.split()) for l in lines[2:] if l.strip()]
        ok, why = verify(kind, points, radius)
        better = Fraction(radius) > Fraction(published[(kind, n)])
        if ok and better and len(points) == n:
            good += 1
        else:
            bad += 1
            print(f"{name}: FAILED ({why}; count {len(points)}; beats published: {better})")
    print(f"{good} packings verified: valid, {''}complete, and larger than the published radius. {bad} failed.")
    return bad == 0


if __name__ == "__main__":
    sys.exit(0 if main(sys.argv[1]) else 1)
