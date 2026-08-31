// constp.cpp -- const(n) = |raw(n)| / prod_{k<2n} k!,  multi-modular.
//
// State = the *avail* bitmask of the reference recursion b(avail,s,d).  s and d are
// determined by avail (s = -sum(avail), d = 2n-1-|avail|), so the memo key is one word, not
// a vector.  Values are held mod a 62-bit prime instead of as bignums; the driver CRTs
// const(n) (not raw(n)) so only ~bits(const)/62 primes are needed.
//   b[0] = 1
//   b[M] = sum_{i in M, ns = s(M)+v_i >= 0}  (-1)^#{j in M : j>i} * g(ns) * b[M ^ 2^i]
//   g(a) = (a+1)(a+2)...(a+n),   v_i = i-(n-1),   raw = b[full]
// Masks are visited in increasing numeric order: every submask is strictly smaller,
// so one flat array suffices and is filled in place.
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <chrono>
using namespace std;
typedef uint64_t u64; typedef unsigned __int128 u128;
int n,m; u64 Q;
static inline u64 mulm(u64 a,u64 b){ return (u64)((u128)a*b % Q); }
static inline u64 addm(u64 a,u64 b){ a+=b; return a>=Q?a-Q:a; }
static inline u64 subm(u64 a,u64 b){ return a>=b?a-b:a+Q-b; }
int main(int argc,char**argv){
    n=atoi(argv[1]); m=2*n-1;
    int NQ=argc-2;
    vector<u64> primes; for(int i=0;i<NQ;i++) primes.push_back(strtoull(argv[2+i],0,10));
    u64 T=1ull<<m;
    int Amax=n*(n-1)/2;
    vector<u64> b(T);
    for(int pi=0;pi<NQ;pi++){
        Q=primes[pi];
        vector<u64> g(Amax+1);
        for(int a=0;a<=Amax;a++){ u64 p=1%Q; for(int t=1;t<=n;t++) p=mulm(p,(u64)(a+t)%Q); g[a]=p; }
        auto t0=chrono::steady_clock::now();
        b[0]=1%Q;
        for(u64 M=1;M<T;M++){
            // s(M) = -sum_{i in M} (i-(n-1)) = (n-1)*popcount(M) - sum of set-bit indices
            int si=0; { u64 Y=M; while(Y){ si+=__builtin_ctzll(Y); Y&=Y-1; } }
            int s=(n-1)*__builtin_popcountll(M)-si;
            u64 t=0; u64 X=M;
            while(X){ u64 lo=X&-X; int i=__builtin_ctzll(lo); X^=lo;
                int ns=s+(i-(n-1)); if(ns<0) continue;
                u64 v=b[M^lo]; if(!v) continue;
                u64 w=mulm(g[ns],v);
                if(__builtin_popcountll(M>>(i+1))&1) t=subm(t,w); else t=addm(t,w);
            }
            b[M]=t;
        }
        double sec=chrono::duration<double>(chrono::steady_clock::now()-t0).count();
        printf("%llu %llu\n",(unsigned long long)Q,(unsigned long long)b[T-1]);
        fflush(stdout);
        fprintf(stderr,"  q=%llu  %.2fs\n",(unsigned long long)Q,sec);
    }
    return 0;
}
