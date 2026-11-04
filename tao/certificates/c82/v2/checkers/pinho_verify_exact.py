"""Independent exact/Arb certificate verifier, written separately from search.

Integrates the positive part using midpoint Taylor quadrature through order 2,
with a rigorous fourth-derivative remainder, on exact dyadic panels. Whole-panel
triangle bounds prove sign. Crossing panels use the positive-part Lipschitz
bound with a rigorous second derivative. No search/upstream imports.
"""
import argparse, hashlib, json, time
from fractions import Fraction
from flint import arb, acb, fmpz_poly, ctx

def rational(x):
 q=Fraction(x);return arb(q.numerator)/q.denominator

def check(path,target,tolerance,precision=192):
 ctx.prec=precision;ctx.threads=1;start=time.time()
 raw=open(path,'rb').read();cert=json.loads(raw)
 den=[];num=[]
 for rec in cert['denominator']:
  c=rec['coeffs'];m=rec['multiplicity']
  assert len(c)>1 and all(type(x) is int for x in c) and abs(c[0])==1
  assert type(m) is int and m>0
  den.append((fmpz_poly(c[::-1]),Fraction(m)))
 for rec in cert['numerator']:
  c=rec['coeffs'];q=Fraction(rec['q_rational'])
  assert len(c)>1 and all(type(x) is int for x in c) and c[0]!=0 and q>=0
  if q:num.append((fmpz_poly(c[::-1]),q))
 D=sum(p.degree()*m for p,m in den);A=sum(p.degree()*q for p,q in num)
 assert 0<=A<D
 for s,_ in num:
  for r,_ in den: assert s.gcd(r).degree()==0
 T=fmpz_poly([0,1,-1]);groups=[];mahler=arb(0);constant=arb(0)
 for sign,group in [(1,num),(-1,den)]:
  for p,q in group:
   cp=fmpz_poly([0])
   for c in reversed(p.coeffs()):cp=cp*T+c
   rr=cp.complex_roots();assert sum(m for r,m in rr)==cp.degree()
   roots=[r for r,m in rr for _ in range(m)]
   w=rational(sign*q/D);lc=arb(abs(int(cp.leading_coefficient()))).log()
   constant+=w*lc;groups.append((w,roots,[abs(r) for r in roots]))
   if sign<0:
    mahler-=w*(lc+sum((abs(r).max(1).log() for r in roots),arb(0)))
 print('exact checks passed; D',D,'A/D',rational(A/D),'Mahler',mahler,flush=True)
 tau=2*arb.pi();tol=rational(tolerance);counts=[0,0,0];splits=0;maxdepth=0
 total=arb(0);stack=[(0,1)];last=time.time()
 while stack:
  k,depth=stack.pop();maxdepth=max(depth,maxdepth)
  width=arb(1)/(1<<depth);h=width/2;mid=arb(2*k+1)/(1<<(depth+1))
  u=acb(0,tau*mid).exp();delta=tau*h
  upper=constant;lower=constant;separated=True;ds=[];possible=True
  for w,rs,ars in groups:
   d=[abs(u-r) for r in rs];ds.append(d)
   lo=arb(1);hi=arb(1);good=True
   for v in d:
    hi*=v+delta
    if v-delta>0:lo*=v-delta
    else:good=False;separated=False
   if w>0:
    upper+=w*hi.log()
    if good:lower+=w*lo.log()
   else:
    if not good:possible=False;break
    upper+=w*lo.log();lower+=w*hi.log()
  if possible and upper<0:counts[0]+=1;continue
  accept=False
  if possible and separated:
   gc=constant;d1=acb(0);d2=acb(0);m2=arb(0);m4=arb(0)
   for (w,rs,ars),distances in zip(groups,ds):
    product=arb(1);first=acb(0);second=acb(0);bound2=arb(0);bound4=arb(0)
    for r,ar,d in zip(rs,ars,distances):
     product*=d;v=u-r;inv=1/v
     first+=inv;second+=r*inv*inv
     low=d-delta;invlow2=1/(low*low)
     bound2+=ar*invlow2
     bound4+=ar*(1+4*ar+ar*ar)*invlow2*invlow2
    gc+=w*product.log();d1+=w*first;d2+=w*second
    m2+=abs(w)*bound2;m4+=abs(w)*bound4
   slope=abs((acb(0,tau)*u*d1).real);curvature=(tau*tau*u*d2).real
   m2*=tau*tau;m4*=tau**4
   reach=slope*h+m2*h*h/2
   if lower>0 or gc-reach>0:
    err=m4*width**5/1920
    value=width*gc+curvature*width**3/24+err;kind=1
   elif gc+reach<0:
    value=arb(0);err=arb(0);kind=0
   else:
    err=slope*width*width/4+m2*width**3/24
    value=width*gc.max(0)+err;kind=2
   if err<tol*width:
    total+=value.upper();counts[kind]+=1;accept=True
  if not accept:
   assert depth<55,('unresolved panel',k,depth)
   stack.append((2*k+1,depth+1));stack.append((2*k,depth+1));splits+=1
  if time.time()-last>60:
   print('progress counts',counts,'splits',splits,'t',float(mid),'seconds',round(time.time()-start),flush=True);last=time.time()
 integral=(2*total).upper();B=mahler+integral;passed=bool(B.upper()<rational(target))
 return dict(certificate_sha256=hashlib.sha256(raw).hexdigest(),target=target,passed=passed,D=int(D),degree_A=str(A),mahler=str(mahler),integral_upper=str(integral),bound_upper=str(B.upper()),counts=counts,splits=splits,max_depth=maxdepth,precision=precision,tolerance=tolerance,seconds=time.time()-start)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('certificate');p.add_argument('--target',default='0.2534350402');p.add_argument('--tol',default='1e-10');p.add_argument('--precision',type=int,default=192);p.add_argument('--out',default='evidence/verification.json');a=p.parse_args()
 result=check(a.certificate,a.target,a.tol,a.precision)
 open(a.out,'w').write(json.dumps(result,indent=2));print(json.dumps(result,indent=2),flush=True)
 raise SystemExit(0 if result['passed'] else 1)
