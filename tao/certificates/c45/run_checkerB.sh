#!/bin/bash
# Checker B (independently written): recompute base histograms with checkerB/histB and evaluate with checkerB/evalB.py.
set -euo pipefail
cd "$(dirname "$0")"
TH=${1:-$(nproc 2>/dev/null || sysctl -n hw.ncpu)}
mkdir -p outB
${CXX:-g++} -O3 -march=native -std=c++17 -pthread checkerB/histB.cpp -o outB/histB
H=outB/histB
$H --threads=$TH 3 5 7 11 13 17 19 29 31 37 41 61 73 > outB/base_griego13.txt
$H --threads=$TH 3 5 7 11 13 17 19 29 31 37 41 61 73 257 > outB/base_yoo14.txt
$H --threads=$TH $(cat specs/primes_final_base.txt) > outB/base_final.txt
for f in base_griego13 base_yoo14 base_final; do cmp outB/$f.txt data/$f.txt && echo "$f: recomputed histogram identical to data/$f.txt"; done
for s in griego yoo final; do echo "== $s"; python3 -I checkerB/evalB.py specs/B_$s.json | tail -3; done
