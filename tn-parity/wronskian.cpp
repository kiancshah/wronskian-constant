// wronskian.cpp -- ranked-automaton layered DP for
//     T_n = sum_{sigma in Phi_n} sgn(sigma) prod_k floor(E_k(sigma)/n)
// State space: subsets S of V_n = {-(n-1)..n-1}, addressed by their rank in the
// regular language supp(f_1) ("alternating sets"), one layer |S|=k at a time.
// Memory = 2 * (largest layer) * sizeof(value)   instead of a hash map over F_{2n+1} keys.
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <string>
#include <chrono>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;
typedef uint64_t u64; typedef int64_t i64; typedef uint32_t u32; typedef uint8_t u8;

static int n, NPOS, MAXK;
static const int SZ[4]  = {0,1,1,2};      // '.' '+' '-' '*'
static const int SGNL[4]= {0,1,-1,0};
static int STEP[4][4];                     // STEP[q][c], q in 1..3 ; -1 = dead
static inline bool acc(int q){ return q==1||q==2; }

// C[j][q][b] = # valid completions of positions j..NPOS-1 from state q using exactly b elements
static vector<vector<vector<u64>>> C;
// CUM[j][q][c][b] = sum over letters c' < c of C[j+1][STEP[q][c']][b-SZ[c']]
static vector<vector<vector<vector<u64>>>> CUM;

static void build_tables(){
    memset(STEP,-1,sizeof STEP);
    STEP[1][0]=1; STEP[1][1]=2; STEP[1][2]=3; STEP[1][3]=2;   // D
    STEP[2][0]=2; STEP[2][1]=-1;STEP[2][2]=3; STEP[2][3]=2;   // U
    STEP[3][0]=3; STEP[3][1]=2; STEP[3][2]=-1;STEP[3][3]=-1;  // N
    C.assign(NPOS+1, vector<vector<u64>>(4, vector<u64>(MAXK+1,0)));
    for(int q=1;q<4;q++) C[NPOS][q][0] = acc(q)?1:0;
    for(int j=NPOS-1;j>=1;j--) for(int q=1;q<4;q++) for(int b=0;b<=MAXK;b++){
        u64 s=0;
        for(int c=0;c<4;c++){ int q2=STEP[q][c]; if(q2<0) continue; if(b<SZ[c]) continue; s+=C[j+1][q2][b-SZ[c]]; }
        C[j][q][b]=s;
    }
    CUM.assign(NPOS+1, vector<vector<vector<u64>>>(4, vector<vector<u64>>(5, vector<u64>(MAXK+1,0))));
    for(int j=1;j<NPOS;j++) for(int q=1;q<4;q++) for(int b=0;b<=MAXK;b++){
        u64 s=0;
        for(int c=0;c<4;c++){ CUM[j][q][c][b]=s; int q2=STEP[q][c];
            if(q2>=0 && b>=SZ[c]) s+=C[j+1][q2][b-SZ[c]]; }
        CUM[j][q][4][b]=s;
    }
    // position NPOS-1..: CUM[NPOS] unused; guard
    for(int q=1;q<4;q++) for(int c=0;c<5;c++) for(int b=0;b<=MAXK;b++) CUM[NPOS][q][c][b]=0;
    // fix: position j=NPOS-1 uses C[NPOS][.][.] which is defined. loop above covers j<NPOS. ok
}
static u64 layer_size(int k){
    u64 s=0; if(k<=MAXK) s+=C[1][1][k];
    if(k>=1) s+=C[1][2][k-1];
    return s;
}

// ---------------- DP ----------------
static u32 MODP;
static int SPLITD=5;                 // modulus
static inline u32 addm(u32 a,u32 b){ a+=b; return a>=MODP?a-MODP:a; }

template<class V>
struct Solver {
    vector<V> prev, cur;
    int k;
    u64 escapes=0;               // targets outside the language that received a nonzero push
    // per-thread scratch
    // recursion state (thread-local via struct)
    struct TL {
        vector<int> cs,qs,used,rem;
        vector<u64> preK,preK1;
        vector<u64> sufR;        // (NPOS+1)*4
        vector<char> sufOK;
        u64 mask; int A;
        vector<u64> srcIdx; vector<int> srcSgn; vector<char> srcOK;
    };
    void leaf(TL&t,u64 idx,const V*P,V*Cu){
        // t.cs[0..NPOS-1], t.qs[0..NPOS], t.used[0..NPOS] filled
        int NP=NPOS;
        for(int i=0;i<=NP;i++) t.rem[i]=k-t.used[i];
        // suffix ranks with ORIGINAL remaining budgets
        for(int q=1;q<4;q++){ t.sufR[NP*4+q]=0; t.sufOK[NP*4+q]=acc(q)?1:0; }
        for(int i=NP-1;i>=1;i--) for(int q=1;q<4;q++){
            int c=t.cs[i]; int q2=STEP[q][c];
            if(q2<0 || !t.sufOK[(i+1)*4+q2]){ t.sufOK[i*4+q]=0; t.sufR[i*4+q]=0; continue; }
            t.sufOK[i*4+q]=1;
            t.sufR[i*4+q]=CUM[i][q][c][t.rem[i]] + t.sufR[(i+1)*4+q2];
        }
        int cnt=0; i64 tot=0;
        u64 mask=t.mask;
        // enumerate removable elements
        for(int j=0;j<NP;j++){
            int c=t.cs[j]; if(c==0) continue;
            int outs[2][2]; int nouts=0;   // {newletter, value}
            if(j==0){ if(c==1){ outs[nouts][0]=0; outs[nouts][1]=0; nouts++; } }
            else {
                if(c==1){ outs[nouts][0]=0; outs[nouts][1]= j; nouts++; }
                else if(c==2){ outs[nouts][0]=0; outs[nouts][1]=-j; nouts++; }
                else { outs[nouts][0]=2; outs[nouts][1]= j; nouts++;
                       outs[nouts][0]=1; outs[nouts][1]=-j; nouts++; }
            }
            for(int e=0;e<nouts;e++){
                int nc=outs[e][0], val=outs[e][1];
                int q=t.qs[j];
                u64 r; bool ok;
                if(j==0){
                    // position 0: letters {0=out(->D,size0), 1=in(->U,size1)}
                    int q2 = (nc==0)?1:2;
                    u64 base = (nc==0)?0:C[1][1][k-1];
                    ok = t.sufOK[1*4+q2];
                    r = ok ? base + t.sufR[1*4+q2] : 0;
                } else {
                    int q2=STEP[q][nc];
                    ok = (q2>=0) && t.sufOK[(j+1)*4+q2];
                    int b = (k-1) - t.used[j];
                    r = ok ? t.preK1[j] + CUM[j][q][nc][b] + t.sufR[(j+1)*4+q2] : 0;
                }
                if(!ok) continue;                      // source not in language => f=0
                int bit = val + (n-1);
                int par = __builtin_popcountll(mask >> (bit+1)) & 1;
                u32 v = P[r];
                if(!v) continue;
                tot += par ? -(i64)v : (i64)v;
                cnt++;
            }
        }
        i64 md=(i64)MODP;
        i64 g = 1 + (i64)(t.A / n);
        i64 res = ((tot % md) * (g % md)) % md; if(res<0) res+=md;
        Cu[idx]=(V)res;
    }
    void rec(TL&t,int pos,u64 idxK,u64 idxK1,const V*P,V*Cu){
        int NP=NPOS;
        if(pos==NP){ if(acc(t.qs[NP]) && t.used[NP]==k) leaf(t,idxK,P,Cu); return; }
        if(pos==0){
            for(int c=0;c<2;c++){
                int q2=(c==0)?1:2, sz=c;
                if(t.used[0]+sz>k) continue;
                t.cs[0]=c; t.qs[1]=q2; t.used[1]=t.used[0]+sz;
                t.preK[1]=idxK + (c==1? C[1][1][k] : 0);
                t.preK1[1]=idxK1 + (c==1? (k>=1?C[1][1][k-1]:0) : 0);
                u64 om=t.mask; int oA=t.A;
                if(c==1) t.mask |= 1ull<<(n-1);
                rec(t,1,t.preK[1],t.preK1[1],P,Cu);
                t.mask=om; t.A=oA;
            }
            return;
        }
        int q=t.qs[pos];
        for(int c=0;c<4;c++){
            int q2=STEP[q][c]; if(q2<0) continue;
            int nu=t.used[pos]+SZ[c]; if(nu>k) continue;
            if(C[pos+1][q2][k-nu]==0) continue;
            t.cs[pos]=c; t.qs[pos+1]=q2; t.used[pos+1]=nu;
            u64 a = idxK  + CUM[pos][q][c][k-t.used[pos]];
            u64 b = (k>=1)? idxK1 + CUM[pos][q][c][(k-1)-t.used[pos]] : 0;
            t.preK[pos+1]=a; t.preK1[pos+1]=b;
            u64 om=t.mask; int oA=t.A;
            if(c==1||c==3) t.mask |= 1ull<<((n-1)+pos);
            if(c==2||c==3) t.mask |= 1ull<<((n-1)-pos);
            t.A += SGNL[c]*pos;
            rec(t,pos+1,a,b,P,Cu);
            t.mask=om; t.A=oA;
        }
    }
    void init_tl(TL&t){
        t.cs.assign(NPOS+1,0); t.qs.assign(NPOS+1,0); t.used.assign(NPOS+1,0); t.rem.assign(NPOS+1,0);
        t.preK.assign(NPOS+1,0); t.preK1.assign(NPOS+1,0);
        t.sufR.assign((NPOS+1)*4,0); t.sufOK.assign((NPOS+1)*4,0);
        t.mask=0; t.A=0;
    }
    // collect DFS prefixes at depth `split` so layers can be filled in parallel
    void collect(TL t,int pos,u64 a,u64 b,int split,vector<TL>&outT,vector<pair<u64,u64>>&outI,vector<int>&outP){
        if(pos==split || pos==NPOS){ outT.push_back(t); outI.push_back({a,b}); outP.push_back(pos); return; }
        if(pos==0){
            for(int c=0;c<2;c++){ int q2=(c==0)?1:2, sz=c;
                if(t.used[0]+sz>k) continue;
                TL u=t; u.cs[0]=c; u.qs[1]=q2; u.used[1]=t.used[0]+sz;
                u64 na=a+(c==1? C[1][1][k]:0), nb=(k>=1)? b+(c==1? C[1][1][k-1]:0) : 0;
                u.preK[1]=na; u.preK1[1]=nb; if(c==1) u.mask|=1ull<<(n-1);
                collect(u,1,na,nb,split,outT,outI,outP); }
            return;
        }
        int q=t.qs[pos];
        for(int c=0;c<4;c++){ int q2=STEP[q][c]; if(q2<0) continue;
            int nu=t.used[pos]+SZ[c]; if(nu>k) continue;
            if(C[pos+1][q2][k-nu]==0) continue;
            TL u=t; u.cs[pos]=c; u.qs[pos+1]=q2; u.used[pos+1]=nu;
            u64 na=a+CUM[pos][q][c][k-t.used[pos]];
            u64 nb=(k>=1)? b+CUM[pos][q][c][(k-1)-t.used[pos]] : 0;
            u.preK[pos+1]=na; u.preK1[pos+1]=nb;
            if(c==1||c==3){ u.mask|=1ull<<((n-1)+pos); u.A+=pos; }
            if(c==2||c==3){ u.mask|=1ull<<((n-1)-pos); u.A-=pos; }
            collect(u,pos+1,na,nb,split,outT,outI,outP); }
    }
    u32 run(bool verbose){
        u64 L0=layer_size(0);
        prev.assign(L0,0); prev[0]=1%MODP;   // empty set
        for(k=1;k<=MAXK;k++){
            u64 Lk=layer_size(k);
            cur.assign(Lk,0);
            TL t0; init_tl(t0);
            t0.qs[0]=0; t0.used[0]=0; t0.preK[0]=0; t0.preK1[0]=0; t0.mask=0; t0.A=0;
            int split = SPLITD;
            vector<TL> PT; vector<pair<u64,u64>> PI; vector<int> PP;
            collect(t0,0,0,0,split,PT,PI,PP);
            const V*P=prev.data(); V*Cu=cur.data();
            #pragma omp parallel for schedule(dynamic,1)
            for(size_t i=0;i<PT.size();i++){
                TL t=PT[i];
                if(PP[i]==NPOS){ if(acc(t.qs[NPOS])&&t.used[NPOS]==k) leaf(t,PI[i].first,P,Cu); }
                else rec(t,PP[i],PI[i].first,PI[i].second,P,Cu);
            }
            if(verbose) fprintf(stderr,"  layer %2d: %12llu states\n",k,(unsigned long long)Lk);
            prev.swap(cur);
        }
        return prev.empty()?0:(u32)prev[0];
    }
};
int main(int argc,char**argv){
    if(argc<2){ fprintf(stderr,"usage: %s n [modulus] [-v]\n",argv[0]); return 1; }
    n=atoi(argv[1]);
    MODP = (argc>2)? (u32)atoll(argv[2]) : (u32)n;
    bool verbose = (argc>3 && string(argv[3])=="-v");
    if(const char*e=getenv("SPLIT")) SPLITD=atoi(e);
    NPOS=n; MAXK=2*n-1;
    build_tables();
    u64 tot=0,mx=0; for(int k=0;k<=MAXK;k++){ u64 L=layer_size(k); tot+=L; if(L>mx)mx=L; }
    fprintf(stderr,"n=%d  |supp|=%llu  max layer=%llu  mod %u\n",n,(unsigned long long)tot,(unsigned long long)mx,MODP);
    auto t0=chrono::steady_clock::now();
    u32 r;
    if(MODP<256){ Solver<u8> S; r=S.run(verbose); }
    else        { Solver<u32> S; r=S.run(verbose); }
    auto t1=chrono::steady_clock::now();
    double sec=chrono::duration<double>(t1-t0).count();
    printf("n=%d  T_n mod %u = %u   [%.3f s, peak ~%.3f GB]\n",n,MODP,r,sec, 2.0*mx*((MODP<256)?1:4)/1e9);
    return 0;
}
