import json,os,sys,time,subprocess,math,random,itertools
from fractions import Fraction
from collections import Counter
T0=time.time()
CSRC=r'''
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
static int D,N,KIND; static double HALF;
static int *head=0,*nxt=0; static int NC,NCT,hcap=0,NCAP; static double CS;
static double DEADLINE;
static double now(){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+1e-9*t.tv_nsec;}
static unsigned long long rs=88172645463325252ULL;
static double rnd(){rs^=rs<<13;rs^=rs>>7;rs^=rs<<17;return (rs>>11)*(1.0/9007199254740992.0);}
static double gauss(){double u=rnd()+1e-300,v=rnd();return sqrt(-2*log(u))*cos(6.283185307179586*v);}
static inline int cellc(double v){int c=(int)((v+HALF)/CS); if(c<0)c=0; if(c>=NC)c=NC-1; return c;}
static void build_grid(const double*x,double cs){
  int nc=(int)floor(2*HALF/cs); if(nc<1)nc=1; if(nc>NCAP)nc=NCAP;
  NC=nc; CS=2*HALF/nc; NCT=(D==2)?nc*nc:nc*nc*nc;
  if(NCT>hcap){head=realloc(head,sizeof(int)*NCT);hcap=NCT;}
  for(int c=0;c<NCT;c++)head[c]=-1;
  for(int i=0;i<N;i++){int c=cellc(x[i*D]); c=c*NC+cellc(x[i*D+1]); if(D==3)c=c*NC+cellc(x[i*D+2]); nxt[i]=head[c]; head[c]=i;}
}
#define FOR_PAIRS(...) \
 for(int i=0;i<N;i++){ int cx=cellc(x[i*D]),cy=cellc(x[i*D+1]),cz=(D==3)?cellc(x[i*D+2]):0; \
  int zlo=(D==3)?cz-1:0, zhi=(D==3)?cz+1:0; \
  for(int a=cx-1;a<=cx+1;a++){ if(a<0||a>=NC)continue; for(int b=cy-1;b<=cy+1;b++){ if(b<0||b>=NC)continue; \
   for(int c=zlo;c<=zhi;c++){ if(D==3&&(c<0||c>=NC))continue; int cell=(D==3)?(a*NC+b)*NC+c:a*NC+b; \
    for(int j=head[cell];j>=0;j=nxt[j]){ if(j<=i)continue; __VA_ARGS__ } }}}}
static double walld(const double*p){
  if(KIND==0){double w=1e9; for(int k=0;k<D;k++){double v=HALF-fabs(p[k]); if(v<w)w=v;} return w;}
  double s=0; for(int k=0;k<D;k++) s+=p[k]*p[k]; return 1-sqrt(s);
}
static double energy(const double*x,double R,double*g){
  double E=0; if(g) memset(g,0,sizeof(double)*N*D);
  for(int i=0;i<N;i++){
    if(KIND==0){ double lim=HALF-R; for(int k=0;k<D;k++){double v=fabs(x[i*D+k])-lim; if(v>0){E+=v*v; if(g) g[i*D+k]+=2*v*(x[i*D+k]>0?1:-1);}}}
    else { double s=0; for(int k=0;k<D;k++) s+=x[i*D+k]*x[i*D+k]; s=sqrt(s); double v=s-(1-R);
      if(v>0){E+=v*v; if(g&&s>0) for(int k=0;k<D;k++) g[i*D+k]+=2*v*x[i*D+k]/s;}}
  }
  build_grid(x,2*R); double c2=4*R*R;
  FOR_PAIRS( double d2=0; double dd[3]; for(int k=0;k<D;k++){dd[k]=x[i*D+k]-x[j*D+k]; d2+=dd[k]*dd[k];}
    if(d2<c2){ double d=sqrt(d2); double v=2*R-d; E+=v*v; if(g&&d>0){ double f=2*v/d; for(int k=0;k<D;k++){g[i*D+k]-=f*dd[k]; g[j*D+k]+=f*dd[k];} } } )
  return E;
}
static double ract(const double*x,double Rg){
  double w=1e9;
  for(int i=0;i<N;i++){ double v=walld(x+i*D); if(v<w)w=v; }
  double cs=2*Rg*1.25; build_grid(x,cs); double m2=cs*cs;
  FOR_PAIRS( double d2=0; for(int k=0;k<D;k++){double t=x[i*D+k]-x[j*D+k]; d2+=t*t;} if(d2<m2)m2=d2; )
  double r=sqrt(m2)/2; return r<w?r:w;
}
#define M 8
static double *S[M],*Y[M],RHO[M],AL[M]; static double *G,*GN,*XN,*DV,*Q;
static double dot(const double*a,const double*b,int n){double s=0;for(int i=0;i<n;i++)s+=a[i]*b[i];return s;}
static double lbfgs(double*x,double R,int maxit){
  int n=N*D; double E=energy(x,R,G); int k=0,cnt=0,stall=0;
  for(int it=0;it<maxit;it++){
    if(E<1e-36) break;
    if((it&15)==0 && now()>DEADLINE) break;
    memcpy(Q,G,sizeof(double)*n);
    for(int m=0;m<cnt;m++){int s=(k-1-m+M*4)%M; AL[s]=RHO[s]*dot(S[s],Q,n); for(int i=0;i<n;i++)Q[i]-=AL[s]*Y[s][i];}
    double gam=1; if(cnt>0){int s=(k-1+M*4)%M; gam=dot(S[s],Y[s],n)/dot(Y[s],Y[s],n);}
    for(int i=0;i<n;i++)Q[i]*=gam;
    for(int m=cnt-1;m>=0;m--){int s=(k-1-m+M*4)%M; double be=RHO[s]*dot(Y[s],Q,n); for(int i=0;i<n;i++)Q[i]+=S[s][i]*(AL[s]-be);}
    for(int i=0;i<n;i++)DV[i]=-Q[i];
    double gd=dot(G,DV,n); if(!(gd<0)){ for(int i=0;i<n;i++)DV[i]=-G[i]; gd=-dot(G,G,n); cnt=0;}
    double mx=0; for(int i=0;i<n;i++){double a=fabs(DV[i]); if(a>mx)mx=a;}
    double t=1; if(mx*t>0.3*R) t=0.3*R/mx;
    double En=0; int ok=0;
    for(int ls=0;ls<40;ls++){ for(int i=0;i<n;i++)XN[i]=x[i]+t*DV[i]; En=energy(XN,R,GN); if(En<=E+1e-4*t*gd){ok=1;break;} t*=0.5; }
    if(!ok){ if(cnt==0) break; cnt=0; continue; }
    int s=k%M; for(int i=0;i<n;i++){S[s][i]=XN[i]-x[i]; Y[s][i]=GN[i]-G[i];}
    double sy=dot(S[s],Y[s],n);
    if(sy>1e-300){RHO[s]=1/sy; k++; if(cnt<M)cnt++;} else cnt=0;
    memcpy(x,XN,sizeof(double)*n); memcpy(G,GN,sizeof(double)*n);
    if(E-En<1e-13*E){ if(++stall>30)break;} else stall=0;
    E=En;
  }
  return E;
}
static double *Z;
static void localp(double*y,double*ry,double estop,int maxit){
  double eps=1e-3, rcur=*ry; int n=N*D;
  while(eps>estop && now()<DEADLINE){
    memcpy(Z,y,sizeof(double)*n); double R=rcur*(1+eps); lbfgs(Z,R,maxit); double r=ract(Z,R);
    if(r>rcur*(1+1e-16)){ double gain=r/rcur-1; rcur=r; memcpy(y,Z,sizeof(double)*n); eps=fmin(1e-2,fmax(eps*2,gain*4)); }
    else eps*=0.25;
  }
  *ry=rcur;
}
static void local(double*y,double*ry){localp(y,ry,1e-12,2000);}
static int *CNT;
static void contacts(const double*x,double R){
  memset(CNT,0,sizeof(int)*N);
  for(int i=0;i<N;i++) if(walld(x+i*D)<R*(1+1e-6)) CNT[i]++;
  build_grid(x,2*R*1.0001); double c2=4*R*R*(1+2e-6);
  FOR_PAIRS( double d2=0; for(int k=0;k<D;k++){double t=x[i*D+k]-x[j*D+k]; d2+=t*t;} if(d2<c2){CNT[i]++;CNT[j]++;} )
}
static void randpt(double*p){
  for(;;){ double s=0; for(int k=0;k<D;k++){p[k]=(2*rnd()-1)*HALF; s+=p[k]*p[k];} if(KIND==0||s<1) return; }
}
static void reinsert(double*x,int i,double R){
  build_grid(x,2*R); double best=-1,bp[3]={0,0,0},p[3];
  int ns=(D==2)?3000:5000;
  for(int t=0;t<ns;t++){
    randpt(p); double w=walld(p); if(w<=best) continue;
    int cx=cellc(p[0]),cy=cellc(p[1]),cz=(D==3)?cellc(p[2]):0; double dm=CS*CS;
    int zlo=(D==3)?cz-1:0, zhi=(D==3)?cz+1:0;
    for(int a=cx-1;a<=cx+1;a++){ if(a<0||a>=NC)continue; for(int b=cy-1;b<=cy+1;b++){ if(b<0||b>=NC)continue;
      for(int c=zlo;c<=zhi;c++){ if(D==3&&(c<0||c>=NC))continue; int cell=(D==3)?(a*NC+b)*NC+c:a*NC+b;
        for(int j=head[cell];j>=0;j=nxt[j]){ if(j==i)continue; double d2=0; for(int k=0;k<D;k++){double q=p[k]-x[j*D+k]; d2+=q*q;} if(d2<dm)dm=d2; }}}}
    double sc=sqrt(dm)/2; if(w<sc)sc=w;
    if(sc>best){best=sc; for(int k=0;k<D;k++)bp[k]=p[k];}
  }
  for(int k=0;k<D;k++) x[i*D+k]=bp[k];
}
static void save(const char*fn,double*y,double r){
  char tmp[600]; sprintf(tmp,"%s.tmp",fn); FILE*f=fopen(tmp,"w"); fprintf(f,"%.17g\n",r);
  for(int i=0;i<N;i++){for(int k=0;k<D;k++)fprintf(f,"%.17g ",y[i*D+k]); fprintf(f,"\n");}
  fclose(f); rename(tmp,fn);
}
static void pick(int mode,int nd,int kdel,const double*P,const int*CC,char*del,double*key){
  if(mode==0) for(int i=0;i<nd;i++) key[i]=CC[i]+rnd()*0.5;
  else if(mode==1) for(int i=0;i<nd;i++) key[i]=(i>=nd-kdel)?1:0;
  else if(mode==2) for(int i=0;i<nd;i++) key[i]=-CC[i]+rnd()*0.5;
  else if(mode==3){ int c=(int)(rnd()*nd); for(int i=0;i<nd;i++){double d2=0;for(int k=0;k<D;k++){double t=P[i*D+k]-P[c*D+k];d2+=t*t;} key[i]=-d2;} }
  else if(mode==4) for(int i=0;i<nd;i++) key[i]=rnd();
  else for(int i=0;i<nd;i++) key[i]=log(rnd()+1e-300)/((CC[i]+1.0)*(CC[i]+1.0));
  memset(del,0,nd);
  for(int t=0;t<kdel;t++){int b=-1; for(int i=0;i<nd;i++) if(!del[i]&&(b<0||key[i]>key[b]))b=i; del[b]=1;}
}
int main(int argc,char**argv){
  FILE*f=fopen(argv[1],"r"); int n,nd; double rd;
  if(fscanf(f,"%d %d %d %d %lf",&D,&KIND,&n,&nd,&rd)!=5)return 1;
  double*P=malloc(sizeof(double)*nd*D); for(int i=0;i<nd*D;i++) if(fscanf(f,"%lf",&P[i])!=1)return 1; fclose(f);
  DEADLINE=now()+atof(argv[3]); int seed=atoi(argv[4]), strat=atoi(argv[5]);
  rs^=(unsigned long long)(seed+1)*0x9E3779B97F4A7C15ULL; for(int i=0;i<20;i++)rnd();
  HALF=(KIND==0)?0.5:1.0;
  double FEST=pow((double)nd/n,1.0/D)-1; double FMAX=fmax(2e-2,FEST);
  NCAP=(D==2)?(int)sqrt(4.0*nd)+2:(int)cbrt(4.0*nd)+2;
  nxt=malloc(sizeof(int)*nd); CNT=malloc(sizeof(int)*nd); int nn=nd*D;
  G=malloc(8*nn);GN=malloc(8*nn);XN=malloc(8*nn);DV=malloc(8*nn);Q=malloc(8*nn);Z=malloc(8*nn);
  for(int m=0;m<M;m++){S[m]=malloc(8*nn);Y[m]=malloc(8*nn);}
  int kdel=nd-n; char*del=calloc(nd,1); double*key=malloc(8*nd);
  int*CC=calloc(nd,sizeof(int));
  N=nd; contacts(P,rd); memcpy(CC,CNT,sizeof(int)*nd);
  N=n; double*best=malloc(8*nn),*y=malloc(8*nn);
  double rb=-1; int nscan=0,bestt=-1;
  int*perm=malloc(sizeof(int)*nd); for(int i=0;i<nd;i++)perm[i]=i;
  for(int i=nd-1;i>0;i--){int j=(int)(rnd()*(i+1)); int t=perm[i];perm[i]=perm[j];perm[j]=t;}
  if(strat>=5 && kdel>0){
    int first=strat-5; double t0=now(); double scanend=t0+(DEADLINE-t0)*0.6;
    for(int t=0;;t++){
      if(kdel==1 && t>0){ memset(del,0,nd); del[perm[t%nd]]=1; }
      else {
        int mode=first;
        if(t>0){double u=rnd(); mode=(u<0.1)?0:(u<0.35)?3:(u<0.6)?4:5;}
        pick(mode,nd,kdel,P,CC,del,key);
      }
      int m=0; for(int i=0;i<nd;i++) if(!del[i]){for(int k=0;k<D;k++)y[m*D+k]=P[i*D+k]; m++;}
      double r=ract(y,rd); localp(y,&r,1e-8,500); nscan++;
      if(r>rb){rb=r; memcpy(best,y,8*N*D); bestt=t; save(argv[2],best,rb);}
      if(now()>scanend) break;
    }
  } else {
    int mode=(strat==1)?0:(strat==3)?2:(strat==2)?3:(strat==4)?1:4;
    pick(mode,nd,kdel,P,CC,del,key);
    int m=0; for(int i=0;i<nd;i++) if(!del[i]){for(int k=0;k<D;k++)best[m*D+k]=P[i*D+k]; m++;}
    rb=ract(best,rd);
  }
  local(best,&rb); save(argv[2],best,rb);
  long moves=0,acc[4]={0,0,0,0},tried[4]={0,0,0,0},dt[5]={0,0,0,0,0},da[5]={0,0,0,0,0};
  while(now()<DEADLINE){
    memcpy(y,best,8*N*D); double amp=rb*(0.02+0.6*rnd()*rnd()); double u=rnd(); int typ; int bin=0;
    if(u<0.08){ typ=0; for(int i=0;i<N*D;i++) y[i]+=gauss()*amp; }
    else if(u<0.18){ typ=1; int c=(int)(rnd()*N); double rad=rb*(3+8*rnd()); for(int i=0;i<N;i++){double d2=0;for(int k=0;k<D;k++){double t=y[i*D+k]-best[c*D+k];d2+=t*t;} if(d2<rad*rad) for(int k=0;k<D;k++) y[i*D+k]+=gauss()*amp*2;} }
    else if(u<0.30){ typ=2; contacts(y,rb); int kk=1+(int)(rnd()*3);
      for(int q=0;q<kk;q++){ int bi=-1; for(int t=0;t<40;t++){int i=(int)(rnd()*N); if(bi<0||CNT[i]<CNT[bi])bi=i;} reinsert(y,bi,rb); CNT[bi]=99; }
      for(int i=0;i<N*D;i++) y[i]+=gauss()*rb*0.01; }
    else { typ=3; double fl;
      if(rnd()<0.85) fl=log(3e-5)+rnd()*(log(3e-3)-log(3e-5)); else fl=log(3e-3)+rnd()*(log(FMAX)-log(3e-3));
      double fr=exp(fl);
      bin=(int)floor((log10(fr)+4.5)*2); if(bin<0)bin=0; if(bin>4)bin=4; dt[bin]++;
      for(int i=0;i<N*D;i++) y[i]+=gauss()*rb*0.003;
      double R=rb*(1+fr);
      for(int st=0;st<12&&now()<DEADLINE;st++){ lbfgs(y,R,3000); double r=ract(y,R); if(r>rb*(1+1e-12)) break; R=rb+(R-rb)*0.5; }
    }
    tried[typ]++;
    lbfgs(y,rb*(1-1e-4),3000); double ry=ract(y,rb); local(y,&ry); moves++;
    if(ry>rb){rb=ry; memcpy(best,y,8*N*D); save(argv[2],best,rb); acc[typ]++; if(typ==3)da[bin]++;}
  }
  FILE*sf=fopen("work/stat.txt","a"); if(sf){fprintf(sf,"%d %d %ld %ld %ld %ld %ld %ld %ld %ld %ld %d %d %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld\n",KIND*10+D,N,moves,tried[0],tried[1],tried[2],tried[3],acc[0],acc[1],acc[2],acc[3],nscan,bestt,dt[0],dt[1],dt[2],dt[3],dt[4],da[0],da[1],da[2],da[3],da[4]); fclose(sf);}
  return 0;
}
'''
open('pack.c','w').write(CSRC)
r=subprocess.run(['cc','-O3','-march=native','-o','pack','pack.c','-lm'],capture_output=True,text=True)
if r.returncode: sys.stderr.write(r.stderr[:1000])
cases=json.load(open('/data/packing_cases.json'))
for c in cases:
    for k in ('n','donor_n'): c[k]=int(c[k])
    for k in ('listed_radius','donor_radius'): c[k]=float(c[k])
KD={'csq':(2,0),'cci':(2,1),'scu':(3,0),'ssp':(3,1)}
def load(path,D):
    rows=[]
    for line in open(path):
        s=line.replace(',',' ').split()
        if not s or s[0].startswith('#'): continue
        try: v=[float(t.replace('D','e').replace('d','e')) for t in s]
        except: continue
        rows.append(v)
    cnt=Counter(len(v) for v in rows)
    L=cnt.most_common(1)[0][0]
    rows=[v for v in rows if len(v)==L]
    if L>=D+2: pts=[v[1:1+D] for v in rows]
    elif L==D+1:
        firsts=[v[0] for v in rows]
        isidx=all(abs(a-round(a))<1e-12 for a in firsts) and len(set(firsts))==len(firsts) and max(firsts)>=len(rows)-1
        pts=[v[1:] for v in rows] if isidx else [v[:D] for v in rows]
    else: pts=[v[:D] for v in rows]
    return [p for p in pts if len(p)==D]
def donor_path(c):
    p=str(c.get('donor_file',''))
    if p and os.path.exists(p): return p
    q=os.path.join('/data/packomania',c['kind']+'_coords',os.path.basename(p))
    if p and os.path.exists(q): return q
    return os.path.join('/data/packomania',c['kind']+'_coords','%s%d.txt'%(c['kind'],c['donor_n']))
RAD={}
def radtab(kind):
    if kind not in RAD:
        d={}
        try:
            for line in open('/data/packomania/%s_radius.txt'%kind):
                s=line.split()
                if len(s)>=2:
                    try: d[int(float(s[0]))]=float(s[1].replace('D','e'))
                    except: pass
        except: pass
        RAD[kind]=d
    return RAD[kind]
def verify(kind,pts,rg):
    D,K=KD[kind]; half=0.5 if K==0 else 1.0
    if any(len(p)!=D for p in pts): return None
    strs=[['%.17e'%v for v in p] for p in pts]
    fp=[[float(s) for s in p] for p in strs]
    w=1e9
    for p in fp:
        if K==0: w=min(w,half-max(abs(v) for v in p))
        else: w=min(w,1-math.sqrt(sum(v*v for v in p)))
    cs=2*rg*1.3; nc=max(1,int(2*half/cs)); csz=2*half/nc
    grid={}
    for i,p in enumerate(fp):
        key=tuple(min(nc-1,max(0,int((v+half)/csz))) for v in p); grid.setdefault(key,[]).append(i)
    m2=cs*cs; cand=[]
    offs=list(itertools.product((-1,0,1),repeat=D))
    for key,lst in grid.items():
        for o in offs:
            k2=tuple(a+b for a,b in zip(key,o))
            if k2 not in grid: continue
            for i in lst:
                for j in grid[k2]:
                    if j<=i: continue
                    d2=sum((a-b)**2 for a,b in zip(fp[i],fp[j]))
                    if d2<cs*cs: cand.append((d2,i,j))
                    if d2<m2: m2=d2
    r=min(math.sqrt(m2)/2,w)-3e-15
    if r<=0: return None
    for attempt in range(4):
        rs='%.17e'%r; R=Fraction(rs); ok=True
        for p,ps in zip(fp,strs):
            if K==0:
                for v,s in zip(p,ps):
                    if half-abs(v)<r+1e-9 and abs(Fraction(s))>Fraction(1,2)-R: ok=False
            else:
                if 1-math.sqrt(sum(v*v for v in p))<r+1e-9:
                    if sum(Fraction(s)**2 for s in ps)>(1-R)**2: ok=False
            if not ok: break
        if ok:
            lim=(2*r)**2*(1+1e-8)+1e-14; R4=4*R*R
            for d2,i,j in cand:
                if d2<lim:
                    if sum((Fraction(a)-Fraction(b))**2 for a,b in zip(strs[i],strs[j]))<R4: ok=False;break
        if ok: return r,rs,strs
        r-=1e-14*(10**attempt)
    return None
def emit(results):
    out=[]
    for ci,(r,rs,strs,src) in results.items():
        c=cases[ci]; D=KD[c['kind']][0]
        if r>c['listed_radius'] and len(strs)==c['n'] and all(len(p)==D for p in strs):
            out.append({'kind':c['kind'],'n':c['n'],'radius':rs,'points':strs})
    print(json.dumps({'packings':out})); sys.stdout.flush()
results={}
os.makedirs('work',exist_ok=True)
try: os.remove('work/stat.txt')
except: pass
TIGHT={('scu',101),('scu',340),('ssp',918),('ssp',935)}
ROOM={('ssp',917),('scu',210),('scu',211),('scu',212),('scu',297)}
def cw(ci):
    c=cases[ci]; key=(c['kind'],c['n'])
    if c['n']>=2000: return 4.0
    if key in TIGHT: return 0.5
    if key in ROOM: return 2.0
    return 1.0
crosslog=[]
def contact_counts(kind,fp,r):
    D,K=KD[kind]; half=0.5 if K==0 else 1.0
    n=len(fp); cnt=[0]*n; th=(2*r)**2*(1+4e-6)
    for i,p in enumerate(fp):
        if K==0: w=half-max(abs(v) for v in p)
        else: w=1-math.sqrt(sum(v*v for v in p))
        if w<r*(1+2e-6): cnt[i]+=1
    cs=2*r*1.01; nc=max(1,int(2*half/cs)); csz=2*half/nc; grid={}
    for i,p in enumerate(fp):
        grid.setdefault(tuple(min(nc-1,max(0,int((v+half)/csz))) for v in p),[]).append(i)
    offs=list(itertools.product((-1,0,1),repeat=D))
    for key,lst in grid.items():
        for o in offs:
            k2=tuple(a+b for a,b in zip(key,o))
            if k2 not in grid: continue
            for i in lst:
                for j in grid[k2]:
                    if j<=i: continue
                    if sum((a-b)**2 for a,b in zip(fp[i],fp[j]))<th: cnt[i]+=1; cnt[j]+=1
    return cnt
def cross():
    order=sorted(results,key=lambda ci:-cases[ci]['n'])
    for ci in order:
        c=cases[ci]
        if c['n']>=2000: continue
        for cj in order:
            d=cases[cj]
            if d['kind']!=c['kind'] or not (0<d['n']-c['n']<=6): continue
            if cj not in results or ci not in results: continue
            rj=results[cj][0]
            if rj<=results[ci][0]*(1+1e-13): continue
            try:
                fp=[[float(s) for s in p] for p in results[cj][2]]
                cnt=contact_counts(c['kind'],fp,rj)
                idx=sorted(range(len(fp)),key=lambda i:(cnt[i],random.random()))
                dl=set(idx[:d['n']-c['n']])
                pts=[p for i,p in enumerate(fp) if i not in dl]
                v=verify(c['kind'],pts,rj)
                if v and v[0]>results[ci][0]:
                    crosslog.append('%s%d<-%d:%.1e'%(c['kind'][:2],c['n'],d['n'],v[0]/results[ci][0]-1))
                    results[ci]=(v[0],v[1],v[2],'x%d'%d['n'])
            except Exception as e: crosslog.append('err %s'%e)
info=[]; ALTS={}
for ci,c in enumerate(cases):
    kind=c['kind']; D,K=KD[kind]; n=c['n']
    try: pts=load(donor_path(c),D)
    except Exception as e:
        info.append('loadfail %s%d %s'%(kind,n,e)); continue
    nd=len(pts)
    if nd<n: info.append('short %s%d'%(kind,n)); continue
    open('work/in%d.txt'%ci,'w').write('%d %d %d %d %.17g\n'%(D,K,n,nd,c['donor_radius'])+'\n'.join(' '.join('%.17g'%v for v in p) for p in pts)+'\n')
    try:
        v=verify(kind,pts[:n],c['donor_radius'])
        if v: results[ci]=(v[0],v[1],v[2],'base')
    except Exception as e: info.append('verr %s'%e)
    alts=[]
    try:
        rt=radtab(kind)
        cand=[(rt[m],m) for m in range(n+1,n+80) if m in rt and m!=c['donor_n'] and rt[m]>c['listed_radius']*(1+1e-6)]
        cand.sort(reverse=True)
        for rm,m in cand:
            if len(alts)>=2: break
            path=os.path.join('/data/packomania',kind+'_coords','%s%d.txt'%(kind,m))
            if not os.path.exists(path): continue
            try: q=load(path,D)
            except: continue
            if len(q)<=n: continue
            fn='work/in%d_a%d.txt'%(ci,len(alts))
            open(fn,'w').write('%d %d %d %d %.17g\n'%(D,K,n,len(q),rm)+'\n'.join(' '.join('%.17g'%v for v in p) for p in q)+'\n')
            alts.append((fn,m,rm))
    except Exception as e: info.append('alt %s'%e)
    ALTS[ci]=alts
emit(results)
NP=20; nc=len(cases)
def run_jobs(jobs,jend):
    running=[]
    while jobs or running:
        while jobs and len(running)<NP:
            args,t=jobs.pop(0)
            left=jend-time.time(); t=min(t,left-3)
            if t<1: jobs=[];break
            p=subprocess.Popen(['./pack',args[0],args[1],'%.1f'%t,args[2],args[3]])
            running.append((p,time.time()+t+20))
        time.sleep(0.2)
        nr=[]
        for p,dl in running:
            if p.poll() is None:
                if time.time()>dl: p.kill()
                else: nr.append((p,dl))
        running=nr
def merge(ci,fn,src):
    c=cases[ci]
    if not os.path.exists(fn): return
    try:
        L=open(fn).read().split('\n'); r0=float(L[0]); pts=[[float(t) for t in l.split()] for l in L[1:] if l.strip()]
    except: return
    if len(pts)!=c['n']: return
    if ci in results and r0<=results[ci][0]: return
    v=verify(c['kind'],pts,r0)
    if v and (ci not in results or v[0]>results[ci][0]): results[ci]=(v[0],v[1],v[2],src)
def restart_phase(ph,jend,nj=2):
    try: cross()
    except Exception as e: crosslog.append('cerr %s'%e)
    emit(results)
    tot=NP*(jend-time.time())*0.92
    W={ci:cw(ci) for ci in results}; sw=sum(W.values())
    jobs=[]
    for ci in results:
        c=cases[ci]; D,K=KD[c['kind']]; r,rs,strs,src=results[ci]
        open('work/p%d_%d.txt'%(ph,ci),'w').write('%d %d %d %d %s\n'%(D,K,c['n'],c['n'],rs)+'\n'.join(' '.join(p) for p in strs)+'\n')
        t=tot*W[ci]/sw/nj
        for j in range(nj):
            fn='work/q%d_%d_%d.txt'%(ph,ci,j)
            try: os.remove(fn)
            except: pass
            jobs.append((('work/p%d_%d.txt'%(ph,ci),fn,str(ph*1000003+j*104729+ci+777),'0'),t))
    jobs.sort(key=lambda a:-a[1])
    run_jobs(jobs,jend)
    win=0
    for ci in list(results):
        if time.time()-T0>1700: break
        b=results[ci][0]
        for j in range(nj): merge(ci,'work/q%d_%d_%d.txt'%(ph,ci,j),'ph%d'%ph)
        if results[ci][0]>b: win+=1
    emit(results)
    return win
p1={}; p1src={}
wins=[]
jend1=T0+700
njob=6
W1={ci:cw(ci) for ci in range(nc)}; sw1=sum(W1.values())
tot1=NP*(jend1-time.time())*0.92
jobs=[]; jdon={}
for j in range(njob):
    for ci in range(nc):
        if os.path.exists('work/in%d.txt'%ci):
            inf='work/in%d.txt'%ci; dm=cases[ci]['donor_n']
            al=ALTS.get(ci,[])
            if j>=4 and al:
                inf,dm,_=al[(j-4)%len(al)]
            jdon[(ci,j)]=dm
            fn='work/out%d_%d.txt'%(ci,j)
            try: os.remove(fn)
            except: pass
            jobs.append(((inf,fn,str(j*7919+ci+13),str(5+j%5)),tot1*W1[ci]/sw1/njob))
jobs.sort(key=lambda a:-a[1])
run_jobs(jobs,jend1)
for ci in range(nc):
    for j in range(njob):
        merge(ci,'work/out%d_%d.txt'%(ci,j),'j%d(%d)'%(j,jdon.get((ci,j),0)))
emit(results)
for ci in results: p1[ci]=results[ci][0]; p1src[ci]=results[ci][3]
ends=[T0+1100,T0+1480]
for k,e in enumerate(ends):
    wins.append(restart_phase(k+1,e))
if time.time()-T0<1650:
    try: cross()
    except Exception as e: crosslog.append('cerr %s'%e)
emit(results)
acc=[0]*4; tr=[0]*4; mv=0
DB={}
try:
    for line in open('work/stat.txt'):
        v=[int(t) for t in line.split()]
        mv+=v[2]
        for k in range(4): tr[k]+=v[3+k]; acc[k]+=v[7+k]
        if len(v)>=23:
            g='big' if v[1]>=2000 else ('2d' if v[0]%10==2 else '3d')
            d=DB.setdefault(g,[[0]*5,[0]*5,0,0])
            for b in range(5): d[0][b]+=v[13+b]; d[1][b]+=v[18+b]
            d[2]+=v[2]; d[3]+=v[11]
except Exception as e: info.append('stat %s'%e)
stat={}
for ci,(r,rs,strs,src) in results.items():
    c=cases[ci]; k=c['kind']; st=stat.setdefault(k,[0,0,0.0,0.0])
    st[0]+=1
    if r>c['listed_radius']:
        st[1]+=1; st[2]+=r/c['listed_radius']-1; st[3]+=r/c['donor_radius']-1
sys.stderr.write('p1src '+' '.join('%d:%s'%(cases[ci]['n'],p1src.get(ci,'-')) for ci in results)+'\n')
sys.stderr.write('final src '+' '.join('%d:%s'%(cases[ci]['n'],results[ci][3]) for ci in results)+'\n')
sys.stderr.write('cross '+' '.join(crosslog)+'\n')
sys.stderr.write('moves=%d tried=%s acc=%s winsperphase=%s\n'%(mv,tr,acc,wins))
for g,d in DB.items(): sys.stderr.write('defl %s moves%d scans%d bins(3e-5,1e-4,3e-4,1e-3,3e-3+) %s\n'%(g,d[2],d[3],' '.join('%d/%d'%(d[1][b],d[0][b]) for b in range(5))))
for k,st in stat.items(): sys.stderr.write('%s beat%d gl=%.3e gd=%.3e\n'%(k,st[1],st[2],st[3]))
g=['%s%d<%d:%.1e/%.1e'%(cases[ci]['kind'][:2],cases[ci]['n'],cases[ci]['donor_n'],p1.get(ci,0)/cases[ci]['donor_radius']-1,results[ci][0]/cases[ci]['donor_radius']-1) for ci in results]
sys.stderr.write('gd '+' '.join(g)+'\n')
sys.stderr.write('info %s time %.0f\n'%(info[:4],time.time()-T0))
emit(results)
