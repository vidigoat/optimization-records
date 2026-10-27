// Exact base-histogram DP for the C45 obstruction, version 3 (sharded, parallel, compact states).
// Usage: basedp3 [--threads N] [--final K] [--shards S] q1 q2 ... (processing order)
// Output: lines "nu delta(nu)" with delta(nu) = #{a mod M0 : nu(a) = nu},
//         nu(a) = #{r mod T0 : a != 2^r (mod q) for all q in the list}.
// Stage 1 (first n-K primes): DP over states f: Z_L -> counts of admissible r mod T0 by residue mod L,
//   L = lcm of orders of the not-yet-processed primes. States are canonicalized under the affine maps
//   x -> u x + t of Z_L (u a unit): the set of remaining coset choices is invariant under these maps.
// Stage 2 (last K primes): all their uncovered patterns on Z_Lf are enumerated with multiplicities and
//   nu = sum_{x in pattern} f(x) is accumulated for every state f.
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <string>
#include <thread>
#include <mutex>
#include <algorithm>
#include <map>
#include <atomic>
#include <chrono>
using namespace std;
typedef unsigned __int128 u128;

static long ordmod2(long q){ long x=2%q,d=1; while(x!=1){x=x*2%q;d++;} return d; }
static long gcdl(long a,long b){while(b){long t=a%b;a=b;b=t;}return a;}
static long lcml(long a,long b){return a/gcdl(a,b)*b;}
static string u128str(u128 x){ if(!x) return "0"; string s; while(x){ s+=char('0'+(int)(x%10)); x/=10;} reverse(s.begin(),s.end()); return s; }
static double now(){ return chrono::duration<double>(chrono::steady_clock::now().time_since_epoch()).count(); }

template<class V> static inline uint64_t hashv(const V* v,int L){
  uint64_t h=0x9E3779B97F4A7C15ull;
  const unsigned char* b=(const unsigned char*)v; size_t nb=(size_t)L*sizeof(V); size_t i=0;
  for(;i+8<=nb;i+=8){ uint64_t w; memcpy(&w,b+i,8); h^=w; h*=0xff51afd7ed558ccdull; h^=h>>32; }
  for(;i<nb;i++){ h^=b[i]; h*=0xc4ceb9fe1a85ec53ull; h^=h>>29; }
  h^=h>>33; h*=0xff51afd7ed558ccdull; h^=h>>33; return h;
}
template<class V> static void least_rotation(const V* v,int L,V* out){
  int i=0,j=1,k=0;
  while(i<L && j<L && k<L){
    V a=v[(i+k)%L], b=v[(j+k)%L];
    if(a==b){ k++; continue; }
    if(a>b) i+=k+1; else j+=k+1;
    if(i==j) j++;
    k=0;
  }
  int s=min(i,j);
  for(int t=0;t<L;t++) out[t]=v[(s+t)%L];
}
template<class V> static bool lessvec(const V* a,const V* b,int L){ for(int i=0;i<L;i++){ if(a[i]!=b[i]) return a[i]<b[i]; } return false; }

template<class V> struct Shard {
  int L=0; vector<V> keys; vector<u128> cnt; vector<uint32_t> slots; size_t mask=0; size_t n=0; mutex mu;
  static constexpr uint32_t EMPTY=0xffffffffu;
  void init(int L_){ L=L_; keys.clear(); cnt.clear(); n=0; slots.assign(64,EMPTY); mask=63; }
  void grow(){ size_t s=slots.size()*2; slots.assign(s,EMPTY); mask=s-1; for(size_t i=0;i<n;i++){ uint64_t h=hashv(&keys[i*L],L); size_t p=(h>>8)&mask; while(slots[p]!=EMPTY) p=(p+1)&mask; slots[p]=i; } }
  void add(const V* k,uint64_t h,u128 c){
    lock_guard<mutex> g(mu);
    if((n+1)*2>slots.size()) grow();
    size_t p=(h>>8)&mask;
    while(slots[p]!=EMPTY){ uint32_t i=slots[p]; if(memcmp(&keys[(size_t)i*L],k,L*sizeof(V))==0){ cnt[i]+=c; return; } p=(p+1)&mask; }
    if(n>=0xfffffff0u){ fprintf(stderr,"shard overflow\n"); exit(1); }
    slots[p]=n; keys.insert(keys.end(),k,k+L); cnt.push_back(c); n++;
  }
  void release(){ vector<V>().swap(keys); vector<u128>().swap(cnt); vector<uint32_t>().swap(slots); n=0; }
};

template<class V> int run(vector<long> qs,int nthreads,int K,int NS){
  int n=qs.size(); vector<long> ds(n); long T0=1; u128 M0=1;
  for(int i=0;i<n;i++){ ds[i]=ordmod2(qs[i]); T0=lcml(T0,ds[i]); M0*=qs[i]; }
  vector<long> Ls(n+1,1); for(int i=n-1;i>=0;i--) Ls[i]=lcml(Ls[i+1],ds[i]);
  if(K<1) K=1; if(K>n) K=n; int n1=n-K;
  fprintf(stderr,"T0=%ld M0=%s n=%d final=%d valbytes=%zu\n",T0,u128str(M0).c_str(),n,K,sizeof(V));
  double t0=now();
  vector<Shard<V>>* cur=new vector<Shard<V>>(NS);
  for(auto&s:*cur) s.init(Ls[0]);
  { vector<V> v(Ls[0],(V)(T0/Ls[0])); (*cur)[0].add(v.data(),hashv(v.data(),Ls[0]),1); }
  for(int st=0;st<n1;st++){
    long q=qs[st], d=ds[st]; int L=Ls[st], L2=Ls[st+1];
    vector<int> units; units.push_back(1); for(int u=2;u<L2;u++) if(gcdl(u,L2)==1) units.push_back(u);
    vector<Shard<V>>* nxt=new vector<Shard<V>>(NS); for(auto&s:*nxt) s.init(L2);
    atomic<size_t> next(0);
    auto work=[&](){
      vector<uint32_t> acc(L2); vector<V> tmp(L2),tmp2(L2),best(L2),cand(L2);
      auto emit=[&](u128 c){
        for(int x=0;x<L2;x++){ if(acc[x]>(uint32_t)(V)(~(V)0)){ fprintf(stderr,"value overflow\n"); exit(1);} tmp[x]=(V)acc[x]; }
        least_rotation(tmp.data(),L2,best.data());
        for(size_t ui=1;ui<units.size();ui++){ long u=units[ui]; for(int x=0;x<L2;x++) tmp2[u*x%L2]=tmp[x]; least_rotation(tmp2.data(),L2,cand.data()); if(lessvec(cand.data(),best.data(),L2)) best.swap(cand); }
        uint64_t h=hashv(best.data(),L2);
        (*nxt)[h%NS].add(best.data(),h,c);
      };
      while(true){
        size_t si=next.fetch_add(1); if(si>=(size_t)NS) break;
        auto& sh=(*cur)[si];
        for(size_t i=0;i<sh.n;i++){
          const V* v=&sh.keys[i*(size_t)L]; u128 c=sh.cnt[i];
          if(q-d>0){ fill(acc.begin(),acc.end(),0); for(int x=0;x<L;x++) acc[x%L2]+=v[x]; emit(c*(u128)(q-d)); }
          for(int s=0;s<d;s++){
            // skip if the class s mod d is already empty: then result equals "no removal" (but weight 1)
            fill(acc.begin(),acc.end(),0); for(int x=0;x<L;x++) if(x%d!=s) acc[x%L2]+=v[x]; emit(c);
          }
        }
      }
    };
    vector<thread> th; for(int t=0;t<nthreads;t++) th.emplace_back(work); for(auto&t:th) t.join();
    for(auto&s:*cur) s.release(); delete cur; cur=nxt;
    size_t tot=0; for(auto&s:*cur) tot+=s.n;
    fprintf(stderr,"step %d q=%ld d=%ld L=%d->%d states=%zu  t=%.1fs\n",st+1,q,d,L,L2,tot,now()-t0);
  }
  int Lf=Ls[n1];
  // patterns of final primes n1..n-2 on Z_Lf; the last prime (qz,dz) is applied explicitly per pattern.
  long qz=qs[n-1], dz=ds[n-1];
  map<vector<uint8_t>,u128> pat; { vector<uint8_t> all(Lf,0xff); pat[all]=1; }
  for(int st=n1;st<n-1;st++){
    long q=qs[st], d=ds[st]; map<vector<uint8_t>,u128> np;
    for(auto& kv:pat){
      if(q-d>0) np[kv.first]+=kv.second*(u128)(q-d);
      for(int s=0;s<d;s++){ vector<uint8_t> m=kv.first; for(int x=s;x<Lf;x+=d) m[x]=0; np[m]+=kv.second; }
    }
    pat.swap(np);
  }
  size_t P=pat.size(); vector<uint8_t> pm(P*Lf); vector<uint64_t> pw(P);
  { size_t i=0; u128 wsum=0; for(auto&kv:pat){ memcpy(&pm[i*Lf],kv.first.data(),Lf); pw[i]=(uint64_t)kv.second; wsum+=kv.second*(u128)qz; i++; }
    if(wsum>>63){ fprintf(stderr,"final weight overflow\n"); exit(1);} }
  size_t totst=0; for(auto&s:*cur){ totst+=s.n; vector<uint32_t>().swap(s.slots); }
  fprintf(stderr,"final stage: Lf=%d states=%zu patterns=%zu last=%ld(d=%ld)\n",Lf,totst,P,qz,dz);
  vector<vector<u128>> H(nthreads, vector<u128>(T0+1,0));
  atomic<size_t> next(0);
  auto work2=[&](int tid){
    vector<u128>& h=H[tid]; vector<uint64_t> loc(T0+1,0); vector<int> touched; touched.reserve(T0+1);
    vector<uint16_t> f16(Lf), g(Lf); vector<uint32_t> rem(dz);
    while(true){
      size_t si=next.fetch_add(1); if(si>=(size_t)NS) break;
      auto& sh=(*cur)[si];
      for(size_t i=0;i<sh.n;i++){
        const V* v=&sh.keys[i*(size_t)Lf];
        for(int x=0;x<Lf;x++) f16[x]=(uint16_t)v[x];
        for(size_t p=0;p<P;p++){
          const uint8_t* m=&pm[p*Lf]; uint32_t tot=0;
          for(int x=0;x<Lf;x++){ uint16_t y=m[x]? f16[x]:0; g[x]=y; tot+=y; }
          uint64_t w=pw[p];
          if(qz-dz>0){ if(loc[tot]==0) touched.push_back(tot); loc[tot]+=w*(uint64_t)(qz-dz); }
          for(int sz=0;sz<dz;sz++) rem[sz]=0;
          for(int x=0;x<Lf;x++) rem[x%dz]+=g[x];
          for(int sz=0;sz<dz;sz++){ uint32_t s2=tot-rem[sz]; if(loc[s2]==0) touched.push_back(s2); loc[s2]+=w; }
        }
        u128 c=sh.cnt[i];
        for(int s2:touched){ h[s2]+=c*(u128)loc[s2]; loc[s2]=0; }
        touched.clear();
      }
    }
  };
  { vector<thread> th; for(int t=0;t<nthreads;t++) th.emplace_back(work2,t); for(auto&t:th) t.join(); }
  vector<u128> hist(T0+1,0); u128 tot=0; for(auto&h:H) for(long s=0;s<=T0;s++) hist[s]+=h[s];
  for(long s=0;s<=T0;s++) tot+=hist[s];
  fprintf(stderr,"done t=%.1fs\n",now()-t0);
  if(tot!=M0){ fprintf(stderr,"TOTAL MISMATCH\n"); return 1; }
  for(long s=0;s<=T0;s++) if(hist[s]) printf("%ld %s\n",s,u128str(hist[s]).c_str());
  return 0;
}

int main(int argc,char**argv){
  int nthreads=1; int K=3; int NS=1024; vector<long> qs; int force16=0;
  for(int i=1;i<argc;i++){ string a=argv[i];
    if(a=="--threads") nthreads=atoi(argv[++i]);
    else if(a=="--final") K=atoi(argv[++i]);
    else if(a=="--shards") NS=atoi(argv[++i]);
    else if(a=="--u16") force16=1;
    else qs.push_back(atol(a.c_str())); }
  int n=qs.size(); long T0=1; vector<long> ds(n); for(int i=0;i<n;i++){ ds[i]=ordmod2(qs[i]); T0=lcml(T0,ds[i]); }
  vector<long> Ls(n+1,1); for(int i=n-1;i>=0;i--) Ls[i]=lcml(Ls[i+1],ds[i]);
  int Kc=max(1,min(K,n)); long maxval=T0/Ls[n-Kc];
  if(!force16 && maxval<=255) return run<uint8_t>(qs,nthreads,K,NS);
  return run<uint16_t>(qs,nthreads,K,NS);
}
