// constp2.cpp -- const(n), multi-modular, LAYERED by |avail| with colex ranking.
//
// b[M] = sum_{i in M, ns = s(M)+v_i >= 0} (-1)^#{j in M : j>i} * g(ns) * b[M \ {i}]
//   v_i = i-(n-1),  s(M) = -sum_{i in M} v_i,  g(a) = (a+1)...(a+n),  raw = b[full].
// b[M] only ever reads masks of popcount |M|-1, so two layers suffice: peak memory is
// C(2n-1, n) + C(2n-1, n-1) values instead of 2^(2n-1).  Colex rank of a k-mask with
// bit positions p_1<...<p_k is sum_j C(p_j, j), and deleting p_j gives
//   rank_{k-1} = sum_{i<j} C(p_i,i) + sum_{i>j} C(p_i,i-1) = Pre[j] + Suf[j],
// so one O(k) pass per mask makes every one of its k transitions O(1).
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <chrono>
#include <string>
#include <new>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;
typedef uint64_t u64; typedef uint32_t u32; typedef unsigned __int128 u128;
static int n,m; static u64 Q;
static vector<vector<u64>> BC;                 // binomials
static inline u64 mulm(u64 a,u64 b){ return (u64)((u128)a*b % Q); }

template<class V>
static u64 solve(const vector<u64>&g,int Amax){
    vector<V> prev(1), cur;
    prev[0] = (V)(1%Q);
    fprintf(stderr,"  n=%d q=%llu  %zu-byte values  peak ~%.2f GB\n",
            n,(unsigned long long)Q,sizeof(V),
            (double)(BC[m][n]+BC[m][n-1])*sizeof(V)/1073741824.0);
    int P[64];  u64 Pre[64], Suf[64];
    for(int k=1;k<=m;k++){
        u64 Lk=BC[m][k];
        try { cur.assign(Lk,0); }
        catch(const std::bad_alloc&){
            fprintf(stderr,"FATAL: out of memory at layer %d/%d (need %.2f GB for this layer "
                           "on top of %.2f GB already held). Raise --mem, or set NARROW=1 "
                           "for 4-byte values.\n",
                    k,m,(double)Lk*sizeof(V)/1073741824.0,
                    (double)prev.size()*sizeof(V)/1073741824.0);
            exit(3);
        }
        u64 nchunk = 1;
        #ifdef _OPENMP
        nchunk = (u64)omp_get_max_threads()*8; if(nchunk>Lk) nchunk=Lk?Lk:1;
        #endif
        #pragma omp parallel for schedule(dynamic,1) firstprivate(P,Pre,Suf)
        for(u64 ch=0; ch<nchunk; ch++){
            u64 lo=Lk*ch/nchunk, hi=Lk*(ch+1)/nchunk;
            if(lo>=hi) continue;
            // colex-unrank `lo` into bit positions
            { u64 r=lo; for(int j=k;j>=1;j--){ int p=j-1; while(p+1<=m-1 && BC[p+1][j]<=r) p++; P[j-1]=p; r-=BC[p][j]; } }
            for(u64 idx=lo; idx<hi; idx++){
                int s=0; for(int j=0;j<k;j++) s += (n-1)-P[j];        // s(M) = -sum v_i
                u64 a=0; for(int j=0;j<k;j++){ Pre[j]=a; a+=BC[P[j]][j+1]; }
                u64 b2=0; for(int j=k-1;j>=0;j--){ Suf[j]=b2; b2+=BC[P[j]][j]; }
                u64 t=0;
                for(int j=0;j<k;j++){
                    int ns = s + (P[j]-(n-1)); if(ns<0) continue;
                    u64 src = prev[Pre[j]+Suf[j]]; if(!src) continue;
                    u64 w = mulm(g[ns], src);
                    if((k-1-j)&1){ t = t>=w? t-w : t+Q-w; } else { t+=w; if(t>=Q) t-=Q; }
                }
                cur[idx]=(V)t;
                // next mask in colex order
                if(idx+1<hi){ int j=0; while(j<k && P[j]+1==(j+1<k?P[j+1]:m)) j++;
                    P[j]++; for(int i=0;i<j;i++) P[i]=i; }
            }
        }
        prev.swap(cur); cur.clear(); cur.shrink_to_fit();
        if(getenv("PROGRESS")) fprintf(stderr,"    layer %2d/%d done (%llu states)\n",k,m,(unsigned long long)Lk);
    }
    return (u64)prev[0];
}
int main(int argc,char**argv){
    n=atoi(argv[1]); m=2*n-1;
    bool wide = getenv("WIDE")!=nullptr;                 // 8-byte values (61-bit primes)
    BC.assign(m+1, vector<u64>(m+2,0));
    for(int i=0;i<=m;i++){ BC[i][0]=1; for(int j=1;j<=i;j++) BC[i][j]=BC[i-1][j-1]+(j<=i-1?BC[i-1][j]:0); }
    int Amax=n*(n-1)/2;
    if(argc>2 && string(argv[2])=="--report"){
        double w=(double)(BC[m][n]+BC[m][n-1]);
        printf("n=%d  peak_states=%.0f  wide_GB=%.2f  narrow_GB=%.2f  dense_GB=%.2f\n",
               n,w,w*8/1073741824.0,w*4/1073741824.0,(double)(1ull<<m)*8/1073741824.0);
        return 0;
    }
    for(int pi=2;pi<argc;pi++){
        Q=strtoull(argv[pi],0,10);
        if(!wide && Q>=(1ull<<32)){
            fprintf(stderr,"FATAL: prime %llu needs more than 32 bits but values are 4-byte "
                           "(NARROW). Residues would be truncated silently. Either set WIDE=1, "
                           "or regenerate the list with: primes.py N --narrow\n",
                    (unsigned long long)Q);
            return 4;
        }
        if(Q<(u64)(2*n)){ fprintf(stderr,"FATAL: prime must exceed 2n so W is invertible mod q\n"); return 4; }
        vector<u64> g(Amax+1);
        for(int a=0;a<=Amax;a++){ u64 p=1%Q; for(int t=1;t<=n;t++) p=mulm(p,(u64)(a+t)%Q); g[a]=p; }
        auto t0=chrono::steady_clock::now();
        u64 r = wide ? solve<u64>(g,Amax) : solve<u32>(g,Amax);
        double sec=chrono::duration<double>(chrono::steady_clock::now()-t0).count();
        printf("%llu %llu\n",(unsigned long long)Q,(unsigned long long)r); fflush(stdout);
        fprintf(stderr,"  q=%llu  %.2fs\n",(unsigned long long)Q,sec);
    }
    return 0;
}
