// dense.cpp -- assumption-free reference: full 2^(2n-1) subset DP, no support hypothesis.
// Confirms (a) T_n, (b) supp(f_g) is exactly the 4-state-automaton language.
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <cstring>
using namespace std;
typedef int64_t i64; typedef uint64_t u64;
int n,m; 
static int STEP[4][4];
static bool inLang(u64 S){
    memset(STEP,-1,sizeof STEP);
    STEP[1][0]=1;STEP[1][1]=2;STEP[1][2]=3;STEP[1][3]=2;
    STEP[2][0]=2;STEP[2][1]=-1;STEP[2][2]=3;STEP[2][3]=2;
    STEP[3][0]=3;STEP[3][1]=2;STEP[3][2]=-1;STEP[3][3]=-1;
    int q = (S>>(n-1)&1) ? 2 : 1;
    for(int j=1;j<n;j++){
        int p=(S>>((n-1)+j))&1, mn=(S>>((n-1)-j))&1;
        int c = p&&mn?3:(p?1:(mn?2:0));
        q=STEP[q][c]; if(q<0) return false;
    }
    return q==1||q==2;
}
int main(int argc,char**argv){
    n=atoi(argv[1]); m=2*n-1; u64 T=1ull<<m;
    vector<int> A(T); vector<i64> f(T,0);
    for(u64 S=1;S<T;S++){ u64 lo=S&-S; int i=__builtin_ctzll(lo); A[S]=A[S^lo]+(i-(n-1)); }
    f[0]=1; u64 badSupp=0, langNotSupp=0, suppCnt=0, langCnt=0;
    for(u64 S=1;S<T;S++){
        if(A[S]<0) continue;
        i64 tot=0; u64 X=S;
        while(X){ u64 lo=X&-X; int i=__builtin_ctzll(lo); X^=lo;
            i64 v=f[S^lo]; if(!v) continue;
            tot += (__builtin_popcountll(S>>(i+1))&1) ? -v : v; }
        f[S] = tot ? tot*(1+A[S]/n) : 0;
    }
    for(u64 S=0;S<T;S++){ bool s=f[S]!=0, l=inLang(S); if(s)suppCnt++; if(l)langCnt++;
        if(s&&!l) badSupp++; if(l&&!s) langNotSupp++; }
    printf("n=%2d  T_n=%lld  |supp|=%llu |lang|=%llu  supp\\lang=%llu  lang\\supp=%llu  maxval=",
        n,(long long)f[T-1],(unsigned long long)suppCnt,(unsigned long long)langCnt,
        (unsigned long long)badSupp,(unsigned long long)langNotSupp);
    i64 mx=0; for(u64 S=0;S<T;S++) if(llabs(f[S])>mx) mx=llabs(f[S]);
    printf("%lld\n",(long long)mx);
    return 0;
}
