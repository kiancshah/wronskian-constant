// cert2.cpp -- O(n)-memory, embarrassingly parallel certificate that  T_n = 1.
// Same identity as cert.cpp, but every language-membership test is O(1) via
// incremental automaton states + a backward "suffix accepts" table.
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <chrono>
#ifdef _OPENMP
#include <omp.h>
#endif
typedef uint64_t u64; typedef int64_t i64;
static int n,NPOS,M;
static const int SZ[4]={0,1,1,2};
static int STEP[4][4];
static inline bool accq(int q){return q==1||q==2;}
#define MAXN 48

struct Ctx {                       // everything is stack-resident: O(n) memory
    int c[MAXN];                   // letters c[0..NPOS-1]  (c[0] in {0,1})
    int q[MAXN+1];                 // automaton state before position j
    char liv[(MAXN+1)*4];          // liv[j*4+q] : letters c[j..] from state q end accepting
    u64 mask; int A;
};
// fill liv from position `from` down, given c[] and NPOS
static inline void build_liv(Ctx&x,int from){
    for(int q=1;q<4;q++) x.liv[NPOS*4+q]=accq(q);
    for(int j=NPOS-1;j>=from;j--){
        for(int q=1;q<4;q++){ int q2=STEP[q][x.c[j]];
            x.liv[j*4+q] = (q2>=0 && x.liv[(j+1)*4+q2]); }
    }
}
static inline bool live_with(const Ctx&x,int j,int nc){   // word with c[j]:=nc in language?
    if(j==0){ int q2=(nc==0)?1:2; return x.liv[1*4+q2]!=0; }
    if(x.q[j]<0) return false;                 // word already dead before position j
    int q2=STEP[x.q[j]][nc]; return q2>=0 && x.liv[(j+1)*4+q2]!=0;
}
static inline int cval_k(int k){ if(k&1) return 1; return ((k/2)&1)?-1:1; }
// sum_{s in S} eps(s,S)*c(S\{s})  over s whose removal stays in L ; k=|S|
static inline i64 rhs(const Ctx&x,int k){
    i64 t=0; int cm=cval_k(k-1);
    for(int j=0;j<NPOS;j++){
        int c=x.c[j]; if(c==0) continue;
        int nc[2],val[2],nn=0;
        if(j==0){ nc[0]=0; val[0]=0; nn=1; }
        else if(c==1){ nc[0]=0; val[0]= j; nn=1; }
        else if(c==2){ nc[0]=0; val[0]=-j; nn=1; }
        else { nc[0]=2; val[0]= j; nc[1]=1; val[1]=-j; nn=2; }
        for(int e=0;e<nn;e++){
            if(!live_with(x,j,nc[e])) continue;
            int bit=val[e]+(n-1);
            t += (__builtin_popcountll(x.mask>>(bit+1))&1)? -cm : cm;
        }
    }
    return t;
}
static u64 g_fail=0;
struct Local { u64 L=0,NB=0,fail=0; };

// check the identity at a neighbour R = S + {one element}, R not in L
static inline void check_escape(const Ctx&S,int j,int nc,u64 rmask,int A2,int k2,int addedBit,Local&lo){
    Ctx R=S;                                  // O(n) copy
    R.c[j]=nc; R.mask=rmask; R.A=A2;
    if(j==0) R.q[1]=(nc==0)?1:2; else R.q[j+1]=(R.q[j]<0)?-1:STEP[R.q[j]][nc];
    for(int i=j+1;i<NPOS;i++) R.q[i+1] = (R.q[i]<0) ? -1 : STEP[R.q[i]][R.c[i]];
    build_liv(R,0);
    // canonical owner: R must be visited from its lowest-indexed L-predecessor
    for(int t=0;t<M;t++){
        if(!((R.mask>>t)&1)) continue;
        int v=t-(n-1); int lj = v<0?-v:v; int nc2;
        if(lj==0) nc2=0;
        else { int cc=R.c[lj]; nc2 = (cc==3)? (v>0?2:1) : 0; }
        if(live_with(R,lj,nc2)){ if(t!=addedBit) return; break; }
    }
    lo.NB++;
    i64 r = rhs(R,k2) * (i64)(1 + A2/n);
    if(r!=0){ lo.fail++; if(lo.fail<3) fprintf(stderr,"FAIL escape mask=%llx got %lld\n",(unsigned long long)R.mask,(long long)r); }
}
static void leaf(Ctx&x,int k,Local&lo){
    build_liv(x,0);
    lo.L++;
    if(x.mask){
        i64 r = rhs(x,k) * (i64)(1 + x.A/n);
        if(r != cval_k(k)){ lo.fail++; if(lo.fail<4) fprintf(stderr,"FAIL in-L mask=%llx got %lld want %d\n",(unsigned long long)x.mask,(long long)r,cval_k(k)); }
    }
    for(int i=0;i<M;i++){
        if((x.mask>>i)&1) continue;
        int v=i-(n-1); int A2=x.A+v; if(A2<0) continue;
        int lj = v<0?-v:v; int nc;
        if(lj==0) nc=1; else { int cc=x.c[lj]; nc = (cc==0)? (v>0?1:2) : 3; }
        if(live_with(x,lj,nc)) continue;              // stays in L -> handled as its own leaf
        check_escape(x,lj,nc,x.mask|(1ull<<i),A2,k+1,i,lo);
    }
}
static void walk(Ctx&x,int pos,int k,Local&lo){
    if(pos==NPOS){ if(accq(x.q[NPOS])) leaf(x,k,lo); return; }
    if(pos==0){
        x.c[0]=0; x.q[1]=1; walk(x,1,k,lo);
        x.c[0]=1; x.q[1]=2; x.mask|=1ull<<(n-1); walk(x,1,k+1,lo); x.mask^=1ull<<(n-1);
        return;
    }
    int q=x.q[pos];
    for(int c=0;c<4;c++){ int q2=STEP[q][c]; if(q2<0) continue;
        x.c[pos]=c; x.q[pos+1]=q2;
        u64 om=x.mask; int oA=x.A;
        if(c==1||c==3){ x.mask|=1ull<<((n-1)+pos); x.A+=pos; }
        if(c==2||c==3){ x.mask|=1ull<<((n-1)-pos); x.A-=pos; }
        walk(x,pos+1,k+SZ[c],lo);
        x.mask=om; x.A=oA;
    }
}
int main(int argc,char**argv){
    n=atoi(argv[1]); NPOS=n; M=2*n-1;
    int SPLIT = (argc>2)?atoi(argv[2]):4;   // parallel split depth
    memset(STEP,-1,sizeof STEP);
    STEP[1][0]=1;STEP[1][1]=2;STEP[1][2]=3;STEP[1][3]=2;
    STEP[2][0]=2;STEP[2][1]=-1;STEP[2][2]=3;STEP[2][3]=2;
    STEP[3][0]=3;STEP[3][1]=2;STEP[3][2]=-1;STEP[3][3]=-1;
    auto t0=std::chrono::steady_clock::now();
    // enumerate all prefixes of length SPLIT, then run them in parallel
    struct Pref{ Ctx x; int k; };
    static Pref pf[1<<20]; int np=0;
    { Ctx x; memset(&x,0,sizeof x); x.q[0]=0; x.mask=0; x.A=0;
      // iterative prefix expansion
      struct It{Ctx x;int pos,k;}; static It st[1<<20]; int sp=0;
      st[sp++] = {x,0,0};
      while(sp){ It it=st[--sp];
        if(it.pos==SPLIT || it.pos==NPOS){ pf[np++] = {it.x,it.k}; continue; }
        if(it.pos==0){ It a=it; a.x.c[0]=0; a.x.q[1]=1; a.pos=1; st[sp++]=a;
                       It b=it; b.x.c[0]=1; b.x.q[1]=2; b.x.mask|=1ull<<(n-1); b.pos=1; b.k=it.k+1; st[sp++]=b; continue; }
        int q=it.x.q[it.pos];
        for(int c=0;c<4;c++){ int q2=STEP[q][c]; if(q2<0) continue;
            It a=it; a.x.c[it.pos]=c; a.x.q[it.pos+1]=q2;
            if(c==1||c==3){ a.x.mask|=1ull<<((n-1)+it.pos); a.x.A+=it.pos; }
            if(c==2||c==3){ a.x.mask|=1ull<<((n-1)-it.pos); a.x.A-=it.pos; }
            a.pos=it.pos+1; a.k=it.k+SZ[c]; st[sp++]=a; }
      }
    }
    u64 TL=0,TN=0,TF=0;
    #pragma omp parallel for schedule(dynamic,1) reduction(+:TL,TN,TF)
    for(int i=0;i<np;i++){
        Local lo; Ctx x=pf[i].x;
        if(SPLIT>=NPOS){ if(accq(x.q[NPOS])) leaf(x,pf[i].k,lo); }
        else walk(x,SPLIT,pf[i].k,lo);
        TL+=lo.L; TN+=lo.NB; TF+=lo.fail;
    }
    double s=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
    printf("n=%2d  |L|=%llu  escapes=%llu  failures=%llu  ->  T_%d = 1 %s   [%.2f s, %d tasks, O(n) memory]\n",
        n,(unsigned long long)TL,(unsigned long long)TN,(unsigned long long)TF,n,
        TF?"NOT PROVED":"PROVED",s,np);
    return TF?1:0;
}
