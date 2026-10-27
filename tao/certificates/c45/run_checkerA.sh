#!/bin/bash
# Checker A: recompute base histograms with checkerA/basedp3 and evaluate exactly with checkerA/exact_eval_a.py.
# Usage: ./run_checkerA.sh [threads]   (python3 >= 3.8 needed; the final base needs ~60 GB RAM / 64 cores)
set -euo pipefail
cd "$(dirname "$0")"
TH=${1:-$(nproc 2>/dev/null || sysctl -n hw.ncpu)}
mkdir -p outA
${CXX:-g++} -O3 -std=c++17 -march=native -pthread checkerA/basedp3.cpp -o outA/basedp3
B=outA/basedp3
$B --threads $TH --final 3 17 29 61 37 19 73 3 7 5 31 11 13 41 > outA/base_griego13.txt
$B --threads $TH --final 3 257 17 29 37 19 73 3 7 5 31 11 13 41 61 > outA/base_yoo14.txt
$B --threads $TH --shards 4096 --final ${FINALK:-4} $(cat specs/order_final.txt) > outA/base_final.txt
for f in base_griego13 base_yoo14 base_final; do cmp outA/$f.txt data/$f.txt && echo "$f: recomputed histogram identical to data/$f.txt"; done
for s in griego yoo final; do echo "== $s"; python3 checkerA/exact_eval_a.py specs/A_$s.json | tail -3; done
