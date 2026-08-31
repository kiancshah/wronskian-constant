# certificate.jl -- O(n)-memory proof that  T_n = 1,  for julia-wronskian.
#
#   T_n = sum over late-growing orderings sigma of {2..2n} of  sgn(sigma) * prod_k floor(E_k/n).
# In excursion coordinates this is  f(V_n)  for the subset DP
#   f(S) = [A(S) >= 0] * (1 + A(S) div n) * sum_{s in S} eps(s,S) * f(S \ {s}),   f({}) = 1
# with A(S) = sum(S), eps(s,S) = (-1)^#{t in S : t > s}.
#
# supp(f) is the language L of a 4-state automaton reading the "levels" j = 0..n-1,
# where level j records which of +j, -j lie in S:
#     letters   '.' = neither   '+' = +j only   '-' = -j only   '*' = both
#     START --(0 not in S)--> D        START --(0 in S)--> U
#     D: . ->D  + ->U  * ->U  - ->N        U: . ->U  * ->U  - ->N  + ->dead
#     N: . ->N  + ->U  * ->dead - ->dead   accepting: D, U     (|L| = F_{2n+1})
# and on L the value is  c(S) = 1 if |S| odd, (-1)^(|S|/2) if |S| even; 0 off L.
#
# So T_n = 1 follows by induction on |S| from the purely LOCAL identity
#     c(S) == (1 + A(S) div n) * sum_{s in S} eps(s,S) * c(S \ {s})            (*)
# holding at every S with A(S) >= 0.  (*) is vacuous unless S or one of its
# one-element deletions meets L, so it is enough to check S in L and its escapes.
# Nothing is memoised: this walks L with O(n) working memory.

const STEP = let t = fill(-1, 4, 4)          # t[q, c+1], q in 1..3, letters '.','+','-','*'
    t[1,1]=1; t[1,2]=2; t[1,3]=3; t[1,4]=2   # D
    t[2,1]=2; t[2,2]=-1; t[2,3]=3; t[2,4]=2  # U
    t[3,1]=3; t[3,2]=2; t[3,3]=-1; t[3,4]=-1 # N
    t
end
const LSZ = (0,1,1,2)
acc(q) = q == 1 || q == 2
cval(k) = isodd(k) ? 1 : (iseven(k >> 1) ? 1 : -1)

inlang(c::Vector{Int}, n) = begin
    q = c[1] == 0 ? 1 : 2
    for j in 2:n
        q = STEP[q, c[j]+1]
        q < 0 && return false
    end
    acc(q)
end

"""Rebuild the letter vector of S and evaluate the RHS of (*)."""
function rhs(c::Vector{Int}, mask::UInt64, n::Int, k::Int)
    tot = 0; cm = cval(k-1); scratch = copy(c)
    for j in 0:n-1
        cj = c[j+1]; cj == 0 && continue
        outs = j == 0 ? ((0, 0),) :
               cj == 1 ? ((0, j),) : cj == 2 ? ((0, -j),) : ((2, j), (1, -j))
        for (nc, v) in outs
            scratch[j+1] = nc
            if inlang(scratch, n)
                bit = v + (n-1)
                tot += isodd(count_ones(mask >> (bit+1))) ? -cm : cm
            end
            scratch[j+1] = cj
        end
    end
    tot
end

function certify(n::Int; verbose=true)
    c = zeros(Int, n); nL = 0; nEsc = 0; fails = 0
    function visit(mask::UInt64, A::Int, k::Int)
        nL += 1
        if mask != 0 && rhs(c, mask, n, k) * (1 + A ÷ n) != cval(k)
            fails += 1
        end
        for i in 0:2n-2                                   # escapes: S + one element leaves L
            (mask >> i) & 1 == 1 && continue
            v = i - (n-1); A2 = A + v; A2 < 0 && continue
            lj = abs(v)
            nc = lj == 0 ? 1 : (c[lj+1] == 0 ? (v > 0 ? 1 : 2) : 3)
            old = c[lj+1]; c[lj+1] = nc
            if !inlang(c, n)
                # visit each escape once, from its lowest-indexed L-predecessor
                own = -1
                for t in 0:2n-2
                    ((mask | (UInt64(1) << i)) >> t) & 1 == 1 || continue
                    vt = t - (n-1); lt = abs(vt)
                    o2 = c[lt+1]
                    c[lt+1] = lt == 0 ? 0 : (o2 == 3 ? (vt > 0 ? 2 : 1) : 0)
                    ok = inlang(c, n); c[lt+1] = o2
                    if ok; own = t; break; end
                end
                if own == i
                    nEsc += 1
                    rhs(c, mask | (UInt64(1) << i), n, k+1) * (1 + A2 ÷ n) != 0 && (fails += 1)
                end
            end
            c[lj+1] = old
        end
    end
    function walk(pos::Int, q::Int, mask::UInt64, A::Int, k::Int)
        if pos == n
            acc(q) && visit(mask, A, k); return
        end
        if pos == 0
            c[1] = 0; walk(1, 1, mask, A, k)
            c[1] = 1; walk(1, 2, mask | (UInt64(1) << (n-1)), A, k+1)
            return
        end
        for ch in 0:3
            q2 = STEP[q, ch+1]; q2 < 0 && continue
            c[pos+1] = ch
            m = mask; a = A
            if ch == 1 || ch == 3; m |= UInt64(1) << ((n-1)+pos); a += pos; end
            if ch == 2 || ch == 3; m |= UInt64(1) << ((n-1)-pos); a -= pos; end
            walk(pos+1, q2, m, a, k + LSZ[ch+1])
        end
    end
    walk(0, 0, UInt64(0), 0, 0)
    verbose && println("n=$n  |L|=$nL  escapes=$nEsc  failures=$fails  ->  T_$n = 1 ",
                       fails == 0 ? "PROVED" : "NOT PROVED")
    fails == 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    for n in (isempty(ARGS) ? (3:2:15) : (parse(Int, ARGS[1]):parse(Int, ARGS[1])))
        certify(n)
    end
end
