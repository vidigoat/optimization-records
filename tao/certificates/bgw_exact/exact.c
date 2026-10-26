/* Exact fixed-point bounds for the C16 sums, c_k = q / (n (q n + p)), n = a^2 + b^2, mu = p/q.
 * Every term is floored exactly with 128-bit integer division at scale 2^100, so
 *   sum(floor) <= 2^100 * S <= sum(floor) + count.
 * Uses the quadrant symmetry: Z^2 \ {0} is four rotations of {a >= 1, b >= 0}.
 *   ./exact K p q a_lo a_hi    prints  count S0lo S1lo S2lo   (for 1 <= a in [a_lo, a_hi])     */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef unsigned __int128 u128;
static void pr(u128 x){ char s[64]; int i=63; s[i]=0; if(!x){printf(" 0");return;} while(x){ s[--i]='0'+(int)(x%10); x/=10;} printf(" %s",s+i); }
static uint64_t isqrt(uint64_t x){ uint64_t r=(uint64_t)__builtin_sqrt((double)x); while(r*r>x) r--; while((r+1)*(r+1)<=x) r++; return r; }
int main(int argc,char**argv){
  uint64_t K=strtoull(argv[1],0,10), p=strtoull(argv[2],0,10), q=strtoull(argv[3],0,10);
  uint64_t alo=strtoull(argv[4],0,10), ahi=strtoull(argv[5],0,10);
  const u128 ONE=((u128)1)<<100, N0=ONE*q, N1=ONE*q*q;
  u128 s0=0,s1=0,s2=0; uint64_t cnt=0, K2=K*K;
  for(uint64_t a=alo;a<=ahi;a++){ uint64_t bmax=isqrt(K2-a*a);
    for(uint64_t b=0;b<=bmax;b++){ uint64_t n=a*a+b*b; u128 m=(u128)q*n+p;
      s0+=N0/((u128)n*m); s1+=N1/((u128)n*m*m); s2+=N1/(m*m); cnt++; } }
  printf("%llu",(unsigned long long)cnt); pr(s0); pr(s1); pr(s2); printf("\n"); return 0; }
