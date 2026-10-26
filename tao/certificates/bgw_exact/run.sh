#!/bin/sh
# Rigorous C16 >= 2.156225504429 (Tao's table, entry 16a): exact fixed-point sums at K = 40000, mu = 41/5.
# About 10 minutes on 14 cores. Pass K as the first argument to try another cutoff (K = 4000 takes a second).
set -e
K=${1:-40000}
cc -O3 -o exact exact.c
rm -f part_*.txt
seq 0 55 | xargs -P "${JOBS:-8}" -I{} sh -c "i={}; ./exact $K 41 5 \$((i*$K/56+1)) \$(((i+1)*$K/56)) > part_\$i.txt"
cat part_*.txt > sums_K$K.txt && rm -f part_*.txt
python3 combine.py sums_K$K.txt
