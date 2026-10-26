#!/usr/bin/env python3
"""Mutation controls: each mutated certificate must be REJECTED by both checkers."""
import json, subprocess, sys, os, tempfile, copy
from fractions import Fraction as F
here = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable
base = json.load(open(os.path.join(here, 'certificate-c0446.json')))

def mut_shift_mass_down(d):
    # move mass 1e-4 from the 2nd atom to the atom at c^2: lowers every moment -> A_k fails
    w = [F(x) for x in d['weights']]; w[1] -= F(1, 10**4); w[0] += F(1, 10**4)
    d['weights'] = [str(x) for x in w]; return d
def mut_shift_mass_up(d):
    # move mass 1e-4 from the 2nd atom to the 3rd: raises beta-weighted moments -> B_k fails
    w = [F(x) for x in d['weights']]; w[1] -= F(1, 10**4); w[2] += F(1, 10**4)
    d['weights'] = [str(x) for x in w]; return d
def mut_raise_c(d):
    c = F(4462, 10000); d['c'] = str(c); n = [F(x) for x in d['nodes']]; n[0] = c * c
    d['nodes'] = [str(x) for x in n]; return d
def mut_starve_last_atom(d):
    w = [F(x) for x in d['weights']]; dw = w[-1] * F(99, 100); w[-1] -= dw; w[-2] += dw
    d['weights'] = [str(x) for x in w]; return d
def mut_bad_mass(d):
    w = [F(x) for x in d['weights']]; w[0] += F(1, 10**12); d['weights'] = [str(x) for x in w]; return d
def mut_jung_model(d):
    # same measure, but Pinho's (weaker) two-zero + Jung model: must fail at c = 0.446
    d['model'] = 'pinho_jung'; d['linear_slope'] = '1/100000'
    d['second_tail_t0'] = '999/1000'; d['second_tail_t1'] = '9999/10000'; d['second_tail_k'] = 20000
    return d
def mut_unknown_model(d):
    d['model'] = 'flux_halfdiameter_v2'; return d

muts = [mut_shift_mass_down, mut_shift_mass_up, mut_raise_c, mut_starve_last_atom, mut_bad_mass, mut_jung_model, mut_unknown_model]
ok = True
with tempfile.TemporaryDirectory() as td:
    for m in muts:
        p = os.path.join(td, m.__name__ + '.json')
        json.dump(m(copy.deepcopy(base)), open(p, 'w'))
        for chk in ['check_c68.py', 'check2_c68.py']:
            r = subprocess.run([PY, os.path.join(here, chk), p], capture_output=True, text=True)
            rejected = r.returncode != 0
            print('%-24s %-14s %s' % (m.__name__, chk, 'REJECTED' if rejected else '!!! ACCEPTED'))
            ok &= rejected
    # -O must be refused
    r = subprocess.run([PY, '-O', os.path.join(here, 'check_c68.py'), os.path.join(here, 'certificate-c0446.json')], capture_output=True, text=True)
    print('python -O refused:', r.returncode != 0); ok &= r.returncode != 0
print('ALL MUTATIONS REJECTED' if ok else 'MUTATION CONTROL FAILURE')
sys.exit(0 if ok else 1)
