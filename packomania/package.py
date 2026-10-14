"""Build a Packomania submission from our verified packings, in the site's requested .pck format.

For each entry the best packing we hold (optimised if we have one, otherwise the shortcut) is
checked again in exact arithmetic and written as <kind><n>.pck: radius on line 1, author on line 2,
then one line of coordinates per object. Everything goes into one zip with the author's name and
the date in its file name, as the site asks, plus the checker so anyone can re-verify.

    python3 -m packing.package --author "Vidit Patankar"
"""
import argparse
import datetime
import glob
import json
import os
import zipfile

from . import table, verify

ROOT = os.path.dirname(table.ROOT)
OUT = os.path.join(os.path.dirname(ROOT), "results", "packomania_submission")


def best_packings():
    best = {}
    for d in (os.path.join(ROOT, "packing_harvest"), os.path.join(ROOT, "packing_best")):
        for p in glob.glob(os.path.join(d, "*.json")):
            x = json.load(open(p))
            key = (x["kind"], int(x["n"]))
            if key not in best or verify.exact(x["radius"]) > verify.exact(best[key]["radius"]):
                best[key] = x
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--author", required=True)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    date = datetime.date.today().isoformat()
    name = f"{a.author.replace(' ', '_')}_{date}.zip"
    zpath = os.path.join(OUT, name)
    rows, count = [], 0
    with zipfile.ZipFile(zpath, "w", zipfile.ZIP_DEFLATED) as z:
        for (kind, n), x in sorted(best_packings().items()):
            listed = table.radii(kind)[n]
            ok, why = verify.verify(kind, x["points"], x["radius"])
            if not ok or verify.exact(x["radius"]) <= verify.exact(listed) or len(x["points"]) != n:
                print(f"skipping {kind}{n}: {why}")
                continue
            body = [str(x["radius"]), a.author] + ["  ".join(str(c) for c in p) for p in x["points"]]
            z.writestr(f"{kind}{n}.pck", "\n".join(body) + "\n")
            gain = float(verify.exact(x["radius"]) / verify.exact(listed)) - 1
            rows.append(f"{kind}\t{n}\t{listed}\t{x['radius']}\t{gain:.3e}")
            count += 1
        z.writestr("SUMMARY.tsv", "category\tN\tpublished_radius\tnew_radius\trelative_gain\n" + "\n".join(rows) + "\n")
        z.write(os.path.join(os.path.dirname(os.path.abspath(__file__)), "verify.py"), "verify_packings.py")
    print(f"{count} packings written to {zpath}")


if __name__ == "__main__":
    main()
