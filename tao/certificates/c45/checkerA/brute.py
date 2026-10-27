import sys
from math import lcm
from collections import Counter
from itertools import product
def ordm(q):
    x=2%q;d=1
    while x!=1: x=x*2%q; d+=1
    return d
qs=[int(a) for a in sys.argv[1:]]
ds=[ordm(q) for q in qs]; T=1
for d in ds: T=lcm(T,d)
pw=[[pow(2,r,q) for r in range(T)] for q in qs]
h=Counter()
for a in product(*[range(q) for q in qs]):
    nu=sum(1 for r in range(T) if all(a[i]!=pw[i][r] for i in range(len(qs))))
    h[nu]+=1
for k in sorted(h): print(k,h[k])
