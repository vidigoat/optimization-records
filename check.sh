#!/bin/sh
# Re-check every result in this repository. Needs Python 3.10+, a C compiler, and mpmath.
#   ./check.sh          quick: all Tao certificates at full strength except 16a (smaller cutoff), and the codes (about a minute)
#   ./check.sh --full   also 16a at full strength and every Packomania packing (downloads the release zip)
set -e
cd "$(dirname "$0")"
echo "== Tao's table =="
python3 tao/certificates/elementary.py
python3 tao/certificates/rbm.py
python3 tao/certificates/rbm_cube.py --n 300000 --t 0.013
if [ "$1" = "--full" ]; then python3 tao/certificates/bgw.py --K 4000 --mu 41/5; else python3 tao/certificates/bgw.py --K 1200 --mu 42/5; fi
echo "== codetables.de =="
mkdir -p bin && cc -O3 -o bin/mindist codetables/mindist.c
python3 codetables/check_all.py
if [ "$1" = "--full" ]; then
  echo "== Packomania =="
  [ -f packomania/packings.zip ] || curl -sL -o packomania/packings.zip \
    https://github.com/vidigoat/optimization-records/releases/download/v1.0/packings.zip
  python3 packomania/verify_all.py packomania/packings.zip
fi
