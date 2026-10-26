#!/bin/sh
# Verify both new certificates and replay Pinho's 0.432 certificate with two independent checkers.
set -e
PY=${PYTHON:-python3}
for f in certificate-c0446.json certificate-c044605.json pinho-certificate-jung-0432.json; do
  echo "== $f"; $PY check_c68.py $f | grep -E '"(status|c|min_rel_margin_A|min_rel_margin_B)"'
  $PY check2_c68.py $f | grep -E 'status|minA|minB'
done
$PY mutation_tests.py
