# The algorithm

This note explains what the two solvers do and why they fit in the memory they fit in.
It assumes the recursion in the main README.

The starting point is a memoised recursion `b(avail, s, d)` over the set of steps not
yet used, the running height, and the depth. Implemented directly, in Python or Julia,
it hits a memory ceiling around `n = 13`, and pushing to `n = 15` would need roughly
100 GB. Nothing about that ceiling is intrinsic to the recursion. It is entirely a
matter of how each state is stored.

## Where the memory actually goes

Three independent costs, all of them representational.

**1. The memo key is a vector.** Both `s` and `d` are functions of `avail`. The step
values sum to zero, so `s = -sum(avail)`, and `d = 2n-1-|avail|` because each call
consumes exactly one element. That makes the state a bitmask, which fits in one
`uint64` for every `n` up to 32. Storing it as a tuple or a vector instead costs more
than the values it indexes: measured at `n = 12`, 443 MB of keys against 330 MB of
values. A vector key is also an allocation per branch, which shows up in the wall clock
as well as the resident set.

**2. The whole table stays live.** `b[M]` reads only masks of popcount `|M| - 1`, so
holding every layer at once is holding `2^(2n-1)` values where two adjacent layers would
do.

**3. Values are big integers.** `b` returns `raw(n)`, which carries the whole factor
`W = prod_{k<2n} k!` around with it before that factor is divided out at the very end.
At `n = 15` that is 1274 bits of dead weight in every single entry. The average memo
value at `n = 12` is already 644 bits, and the size grows roughly quadratically in `n`.

There is one thing that does not help, and it is worth being explicit about it. The
number of live states is `0.52 * 2^(2n-1)` for every `n` tested. The weight
`(a+1)...(a+n)` is strictly positive, so no cancellation happens and no sparsity appears.
The win has to come from making each state cheap.

## The three changes

**Bitmask key, filled in numeric order.** Every submask of `M` is numerically smaller
than `M`, so iterating `M = 1, 2, 3, ...` reaches every dependency before its dependent.
That removes the hash table and the recursion stack together: a flat array indexed by the
mask, filled in place. This is `constp.cpp`.

**Layer by popcount, address by colex rank.** Since `b[M]` reads only popcount `|M| - 1`,
two consecutive layers suffice. Inside a layer, a state is located by the colex rank of
its mask. For bit positions `p_1 < ... < p_k` that rank is

```
rank_k(M) = sum_j C(p_j, j)
```

and deleting `p_j` gives

```
rank_{k-1}(M \ {p_j}) = sum_{i<j} C(p_i, i) + sum_{i>j} C(p_i, i-1) = Pre[j] + Suf[j]
```

which is a prefix sum plus a suffix sum. One `O(k)` pass over the mask therefore makes
all `k` of its transitions `O(1)`, so the ranking is asymptotically free, and peak memory
becomes `C(2n-1, n) + C(2n-1, n-1)` instead of `2^(2n-1)`. This is `constp2.cpp`.

**Multi-modular arithmetic.** The recursion uses only `+`, `-` and `*`, and reduction mod
a prime is a ring homomorphism, so running the identical recursion in `Z/qZ` gives
`raw(n) mod q`. With a 61-bit prime and `__int128` multiplies that is 8 bytes per state
and no allocator traffic in the inner loop.

Two details matter here. First, reconstruct `const`, not `raw`. For any `q > 2n` every
factor of `W` is a unit mod `q`, so `const = raw * W^{-1} (mod q)`, and since `raw` is
about twice the size of `const` this halves the number of primes. Second, the primes are
completely independent, so the parallelisation is one prime per core, or one prime per
job array task, with no communication.

`NARROW=1` switches to 30-bit primes and 4-byte residues: half the memory for roughly
twice the number of primes. The solver refuses a prime above 32 bits in narrow mode
rather than truncating it silently.

## What this buys

Measured on two cores. "Straightforward" is the entry count times the average value size
plus the key words, measured exactly up to `n = 12` and extrapolated above; it ignores
allocator and hash table overhead, so it is an underestimate.

| n | straightforward | flat array | layered 8-byte | layered 4-byte | wall clock |
|---|---|---|---|---|---|
| 13 | ~6 GB | 0.26 GB | 0.09 GB | 0.045 GB | 62 s |
| 14 | ~25 GB | 1.0 GB | 0.30 GB | 0.15 GB | 199 s |
| 15 | ~100 GB | 4.0 GB | 1.16 GB | 0.58 GB | ~12 min |
| 16 | ~450 GB | 16 GB | 4.48 GB | 2.24 GB | ~50 min |
| 17 | | 64 GB | 17.4 GB | 8.7 GB | ~3 h |
| 18 | | 256 GB | 67.6 GB | 33.8 GB | ~12 h |

The practical effect is that a 16 GB laptop reaches `n = 16` and a 128 GB cluster node
reaches `n = 18`.

`constp.cpp` is kept even though `constp2.cpp` dominates it on memory, for two reasons:
it is simpler to read, and it is an independent implementation to check against.

## A different weight, and a much better structure

The above is about the exact weight. The mod `p` reduction replaces it with
`floor(A_k/n)`, which is constant on the window `{0, ..., n-1}`, and that changes the
picture completely: the reduced sum does cancel, its support is a regular language on
four states, and the values are all in `{0, +1, -1}`. That is a different piece of work
and lives in `tn-parity/`, with its own README.
