"""Check the record codes: rank and exact minimum distance by enumerating every codeword.

    cc -O3 -o mindist codes/mindist.c && python3 codes/check_all.py
"""
import glob, os, subprocess, sys
here = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(os.path.dirname(here), "mindist")
TABLE = {(27, 10): 12, (64, 10): 38, (81, 10): 50, (82, 10): 51}   # codetables.de lower bounds, 15 Oct 2026
for path in sorted(glob.glob(os.path.join(here, "gf5_*.txt"))):
    G = [list(map(int, l.split())) for l in open(path) if l.strip() and not l.startswith("#")]
    k, n = len(G), len(G[0])
    out = subprocess.run([BIN], input=f"5 {n} {k} 0\n" + " ".join(str(x) for r in G for x in r), capture_output=True, text=True)
    d = int(out.stdout.split()[0])
    print(f"[{n},{k},{d}] over GF(5): table lower bound {TABLE[(n, k)]} -> {'NEW RECORD' if d > TABLE[(n, k)] else 'no'}")
