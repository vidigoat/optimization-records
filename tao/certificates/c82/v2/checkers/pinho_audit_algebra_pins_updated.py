"""Fresh review: exact JSON and modular coprimality audit; no submitted-code imports.

A degree-preserving reduction coprime modulo one prime proves coprimality over Q.
Polynomial arithmetic here is elementary Python on ascending modular arrays.
"""
import json, hashlib, math, sys
from fractions import Fraction
from pathlib import Path
EXPECTED='1125ddf42e9a544c6fc566aa2ebdfb59ec7c35177030c5d50767427e52f841d1'
def unique(pairs):
    out={}
    for k,v in pairs:
        assert k not in out, ('duplicate JSON key',k)
        out[k]=v
    return out

def trim(a):
    while len(a)>1 and a[-1]==0: a.pop()
    return a

def gcd_mod(left,right,p):
    a=[c%p for c in left[::-1]]; b=[c%p for c in right[::-1]]
    assert a[-1] and b[-1]
    while any(b):
        r=a[:]; inv=pow(b[-1],-1,p)
        while any(r) and len(r)>=len(b):
            v=r[-1]*inv%p; shift=len(r)-len(b)
            for j,c in enumerate(b):r[j+shift]=(r[j+shift]-v*c)%p
            trim(r)
        a,b=b,r
    return len(a)-1

def main():
    raw=Path(sys.argv[1]).read_bytes(); sha=hashlib.sha256(raw).hexdigest();assert sha==EXPECTED
    c=json.loads(raw,object_pairs_hook=unique)
    degreeQ=0;degreeA=Fraction(0);b=1; den=[];num=[];records=[]
    for role in ['denominator','numerator']:
        for i,row in enumerate(c[role]):
            a=row['coeffs']
            assert type(a) is list and len(a)>=2 and all(type(v) is int for v in a) and a[0]!=0
            rec=dict(role=role,index=i,name=row.get('name'),degree=len(a)-1,lead=a[0])
            if role=='denominator':
                m=row['multiplicity'];assert type(m) is int and m>0 and abs(a[0])==1
                degreeQ+=(len(a)-1)*m;den.append(row);rec['multiplicity']=m
            else:
                assert type(row['q_rational']) is str
                q=Fraction(row['q_rational']);assert q>=0
                degreeA+=(len(a)-1)*q;b=math.lcm(b,q.denominator)
                rec['exponent']=str(q)
                if q>0:num.append(row)
            records.append(rec)
    assert 0<=degreeA<degreeQ and degreeQ>0
    assert degreeQ==999999999960 and degreeA==Fraction(473650268755560921,524288)
    witnesses=[]
    for s in num:
        for r in den:
            for prime in [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101]:
                if s['coeffs'][0]%prime and r['coeffs'][0]%prime and gcd_mod(s['coeffs'],r['coeffs'],prime)==0:
                    witnesses.append(dict(numerator=s['name'],denominator=r['name'],prime=prime));break
            else: raise AssertionError(('no modular coprimality witness',s['name'],r['name']))
    result=dict(sha256=sha,passed=True,D=degreeQ,degreeA=str(degreeA),slack=str(degreeQ-degreeA),degreeA_over_D=str(degreeA/degreeQ),common_exponent_denominator=b,degreeW=2*b*degreeQ,degreeU=int(2*b*degreeA),coprimality_pairs=len(witnesses),records=records,modular_witnesses=witnesses)
    Path(sys.argv[2]).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['records','modular_witnesses']},indent=2),flush=True)
if __name__=='__main__':main()
