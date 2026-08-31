# tn-parity: the mod-p side (`T_n`), separate from `const(n)`

Not needed to compute `const(n)`. This is about the reduced sum that `wronskian.cpp` evaluates:

    T_n = sum_{sigma in Phi_n} sgn(sigma) * prod_{k=1}^{2n-1} floor(E_k(sigma)/n)   =  1

i.e. `raw(p)/p^(2p-1) = (-1)^p T_p (mod p)`, the prime-free form of the residual conjecture.

## Build

```
make          # gcc/clang + OpenMP   (macOS: make CXX=g++-14, brew install gcc)
make check    # ~1 min: dense reference, certificate and DP all agree
```

| binary | what | memory |
|---|---|---|
| `./wronskian n [mod]` | layered DP over ranked support states | `2*maxlayer(n)` bytes |
| `./cert n [split]` | O(n)-memory certificate that `T_n = 1` | a few kB |
| `./dense n` | full `2^(2n-1)` subset DP, no assumptions | `2^(2n-1)*12` bytes |
| `certificate.jl` | the certificate, in Julia | a few kB |

## What turned up

Unlike the exact `const(n)` weight, this weight is constant on the window `{0..n-1}`, so it
*does* cancel, and `supp(f)` is not merely Fibonacci-sized, it is a **regular language** on
four states reading the levels `j = 0..n-1` (level `j` records which of `+j, -j` lie in `S`):

```
letters:  '.' neither   '+' only +j   '-' only -j   '*' both
START --(0 not in S)--> D        START --(0 in S)--> U
D:  . -> D    + -> U    * -> U    - -> N
U:  . -> U    * -> U    - -> N    + -> DEAD
N:  . -> N    + -> U    * -> DEAD - -> DEAD
accepting: D, U      |L| = F_{2n+1}
```

Verified letter-for-letter against `dense` for n = 2..13: `supp(f) == L` exactly, every
value lies in `{0, +-1}`, and on `L` it equals `1` for `|S|` odd and `(-1)^(|S|/2)` for
`|S|` even.

Two consequences. **Ranking**: a regular language can be ranked, so states become
consecutive integers and a flat byte array replaces the hash map, n=23 goes from ~150 GB
to 0.71 GB. **A local certificate**: since the DP's answer is closed-form, `T_n = 1` follows
by induction from `c(S) == (1 + A(S) div n) * sum_{s in S} eps(s,S) c(S\{s})` at every `S`
with `A(S) >= 0`, which is memoryless, walk `L`, check each state and its escapes, store
nothing. That proves the integer identity rather than a residue, and a failure would be an
explicit counterexample to Obs. 2.

| n | `wronskian` | peak RSS | `cert` (2 threads) | previously |
|---|---|---|---|---|
| 17 | 2.4 s | 3 MB | 11.9 s | 13 s |
| 19 | 17.5 s | 17 MB | 92 s | ~5 min, 7 GB |
| 21 | 128 s | 109 MB | ~13 min est | out of reach |
| 23 | 1014 s | 0.71 GB | ~2 h / cores est | 128 GB node |
