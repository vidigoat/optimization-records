"""mkspec.py name 'B:..' 'A:..' ... -> name.json (spec with clusters of the 54-prime set) and prints the cut C."""
import json, sys
sys.path.insert(0, '.')
from c45eval import ord2, lam_log2, cluster_law
from fractions import Fraction
from math import lcm, prod
name, groups = sys.argv[1], sys.argv[2:]
base = [int(q) for g in groups for q in g[2:].split(',')]
cl = json.load(open('spec54.json'))['clusters']
json.dump({'base': base, 'clusters': cl}, open(name + '.json', 'w'))
cl = {int(k): v for k, v in cl.items()}
M0 = prod(base); T0 = lcm(*[ord2(q) for q in base]); Mc = prod(q for v in cl.values() for q in v)
T = T0 * prod(cl); phi = prod(q - 1 for q in base) * prod(q - 1 for v in cl.values() for q in v)
kappa = Fraction(M0 * Mc, T * phi) / lam_log2()
Hmin = prod(min(cluster_law(d, qs)) for d, qs in cl.items())
print((Fraction(1, 2) / (Hmin * kappa)).__ceil__() - 1)
