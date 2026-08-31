# Why the fast solvers compute the same number

The solvers in `constp.cpp` and `constp2.cpp` look nothing like the recursion in
`reference.py`, so this note argues, step by step, that they compute the same thing.
Nothing in the recursion changed. What changed is the memo key, the order the same table
is filled in, and the ring the arithmetic happens in.

Write `V = {-(n-1), ..., n-1}` for the step values, obtained from `e in {2, ..., 2n}` by
`v = e - n - 1`. Note `|V| = 2n-1` and `sum(V) = 0`.

## 1. `s` is a function of `avail`

At every node, `s` is the sum of the step values already consumed. The consumed set is
`V \ avail`, so

```
s  =  sum(V) - sum(avail)  =  -sum(avail).
```

*Induction.* At the root, `avail = V` and `s = 0 = -sum(V)`. If it holds at `(avail, s)`
and we take `e` with value `v`, then `avail' = avail \ {e}` and
`s' = s + v = -sum(avail) + v = -sum(avail')`.

## 2. `d` is a function of `avail`

`d = (2n-1) - |avail|`, again by induction: `d` starts at 0 with `|avail| = 2n-1`, and
each call increments `d` and removes one element. So the base case `d == N-1` is exactly
`avail == {}`.

Therefore `memo[(avail, s, d)]` and `memo[avail]` are the same table. Both `s` and `d`
are derived data. The key becomes the bitmask `M` of `avail`, with bit `i` corresponding
to `e = i + 2`.

## 3. The sign is the same expression

With `avail` sorted ascending and `e = avail[i]`,

```
len(avail) - 1 - i  =  #{ f in avail : f > e },
```

which in the bitmask is `popcount(M >> (i+1))`. This is identical, not equivalent up to
something.

## 4. The weight is the same expression

`prod(range(ns+1, ns+n+1)) = (ns+1)...(ns+n) =: g(ns)`, tabulated once for
`ns in [0, n(n-1)/2]`, the largest reachable height being the sum of the positive step
values.

So the implemented recursion is, letter for letter, the reference one:

```
b[{}] = 1
b[M]  = sum_{i in M,  ns := s(M) + v_i >= 0}  (-1)^popcount(M >> (i+1)) * g(ns) * b[M \ {i}]
const(n) = |b[V]| / prod_{k=0}^{2n-1} k!
```

## 5. Loop instead of recursion

`M \ {i}` is numerically smaller than `M`, so filling `M = 1, 2, 3, ...` in increasing
integer order reaches every dependency before its dependent. Same values, same partial
order, no recursion stack and no hash table.

## 6. Two layers instead of the whole table

`b[M]` reads only masks of popcount `|M| - 1`. Sweeping by popcount and holding two
consecutive layers therefore computes the same numbers. Inside a layer, a state is
located by its colex rank instead of by its mask, and colex rank is a bijection from
masks of popcount `k` onto `{0, ..., C(2n-1,k)-1}`, so this is relabelling and not a
change of content.

## 7. Arithmetic mod q, then CRT

The recursion uses only `+`, `-` and `*`, and the final division is by the known constant
`W = prod_{k<2n} k!`. Reduction mod a prime `q` is a ring homomorphism, so running the
identical recursion in `Z/qZ` yields `raw(n) mod q`. For `q > 2n` every factor of `W` is
a unit, hence

```
const(n)  ==  raw(n) * W^{-1}   (mod q).
```

Do this for primes `q_1, ..., q_r` and CRT to get `const(n) mod (q_1 ... q_r)`. If
`q_1 ... q_r > 2|const(n)|`, the balanced representative in `(-M/2, M/2]` is `const(n)`
exactly.

The bound used to size the prime list is proved, not estimated:

```
|raw(n)| <= (2n-1)! * g(A_max)^(2n-1),      A_max = n(n-1)/2,
|const(n)| = |raw(n)| / W,
```

since there are at most `(2n-1)!` orderings and each contributes a product of `2n-1`
weights, each at most `g(A_max)`. `cluster/primes.py:bound_bits_const` implements a
sharper version of the same idea in exact integer arithmetic, using only `delta_k >= 0`
and the concavity of `E -> log(E^{falling n})` via Jensen. Both `driver.py` and the
cluster path size their prime lists from that function, so there is no probabilistic
early-termination step anywhere in this repository.

## 8. What was actually checked

- `verify.py` runs `reference.py` and the modular pipeline side by side and asserts
  equality. The reference implementation is the limiting factor here, not the solvers.
- `constp.cpp` (flat array over masks) and `constp2.cpp` (layered, colex-ranked) were
  written independently of each other and agree residue for residue at every `n` and
  every prime tested, and at 1, 2, 8 and 16 threads.
- The values for `n = 1..13` agree with the b-file of OEIS A392714.
- `const(14)` comes out at 241 digits, matching the value in the paper.
- `identity_check.py` expands the defining operator identity symbolically and recovers
  `const(1) = 1` and `const(2) = 2` without using the recursion at all.
