/* Exact C1c score of integer step weights: K * min_{0<=k<=K} A_k / S^2 with A_k = sum_i w_i w_{i+k}.
   Prints the minimising lag, A_min and S as decimal integers (128-bit), then the score as a double.   */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef __int128 i128;
static void pr(i128 x){ char s[64]; int i=63; s[i]=0; if(!x){printf("0");return;} while(x){ s[--i]='0'+(int)(x%10); x/=10;} printf("%s",s+i); }
int main(int argc,char**argv){
  FILE*f=fopen(argv[1],"r"); fseek(f,0,2); long L=ftell(f); rewind(f); char*b=malloc(L+1); fread(b,1,L,f); b[L]=0;
  long K=atol(strstr(b,"\"K\"")+4); char*p=strchr(strstr(b,"weights_int"),'[')+1;
  long n=0,cap=1<<20; long long*w=malloc(cap*8);
  while(*p && *p!=']'){ w[n++]=strtoll(p,&p,10); while(*p==','||*p==' '||*p=='\n') p++; }
  i128 S=0; for(long i=0;i<n;i++) S+=w[i];
  i128 best=-1; long bk=-1;
  #pragma omp parallel for schedule(dynamic,64)
  for(long k=0;k<=K;k++){ i128 A=0; for(long i=0;i+k<n;i++) A+=(i128)w[i]*w[i+k];
    #pragma omp critical
    { if(best<0||A<best||(A==best&&k<bk)){best=A;bk=k;} } }
  printf("K=%ld n=%ld argmin_k=%ld A_min=",K,n,bk); pr(best); printf(" S="); pr(S);
  printf(" score~%.15Lf\n",(long double)K*(long double)best/((long double)S*(long double)S)); }
