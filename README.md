# wronskian-constant

[![build and test](https://github.com/kiancshah/wronskian-constant/actions/workflows/ci.yml/badge.svg)](https://github.com/kiancshah/wronskian-constant/actions/workflows/ci.yml)

Exact computation of `const(n)`, the integer that appears when an alternating
composition of weighted differential operators collapses to a Wronskian.

For weights `w_1, ..., w_{2n}` and `D = d/dx`,

```
sum_{sigma in S_{2n}} sign(sigma) * w_{sigma(1)}D^n o w_{sigma(2)}D^n o ... o w_{sigma(2n)}D^n
        =  const(n) * Wronskian(w_1, ..., w_{2n}) * D^n
```

Everything on the left is an operator of order up to `2n^2`. Almost all of it cancels.
What survives is a single term of order `n`, and its coefficient is a fixed integer
multiple of the Wronskian, independent of the weights. That integer is `const(n)`, and
this repository computes it.

The sequence is [OEIS A392714](https://oeis.org/A392714). The background is
[arXiv:2605.11137](https://arxiv.org/abs/2605.11137).

## A worked example

Take `n = 1`, so two weights and first-order operators. There are two permutations:

```
w1 D o w2 D  =  w1 w2' D  +  w1 w2 D^2
w2 D o w1 D  =  w2 w1' D  +  w1 w2 D^2
```

The second-order parts are identical, so they cancel in the signed sum, and

```
w1 D o w2 D  -  w2 D o w1 D  =  (w1 w2' - w2 w1') D  =  Wronskian(w1, w2) * D
```

which is `const(1) = 1`.

At `n = 2` the same thing happens with 24 terms and operators of order up to 8. Every
order except `D^2` cancels, and the coefficient of `D^2` comes out at exactly twice the
Wronskian, so `const(2) = 2`. Both are checked symbolically by `identity_check.py`,
which expands the operators rather than using the recursion:

```
$ python3 identity_check.py 2
n = 1
  orders present besides D^1: none
  coefficient of D^1, divided by the Wronskian: 1
n = 2
  orders present besides D^2: none
  coefficient of D^2, divided by the Wronskian: 2
```

That approach costs `(2n)!` operator compositions and stops being usable almost
immediately. The rest of this repository is about getting further.

## What the numbers look like

`const(n)` grows superexponentially. The first terms are `1, 2, 90, 586656,
1915103977500`, and after that they are easier to describe by size:

| n | digits | n | digits |
|---|---|---|---|
| 5 | 13 | 11 | 127 |
| 6 | 22 | 12 | 160 |
| 7 | 35 | 13 | 198 |
| 8 | 52 | 14 | 241 |
| 9 | 73 | 15 | 289 |
| 10 | 98 | | |

Exact values are in [`data/const.txt`](data/const.txt), one per line.

## Quick start

Needs a C++17 compiler and Python 3.8 or later. No libraries, no build system beyond
`make`.

```
make                            # builds constp and constp2
python3 driver.py 12            # print const(1) through const(12), a few seconds
python3 verify.py               # run the test suite
```

Larger `n` is the same command with more patience:

```
python3 driver.py 15            # about 12 minutes, 1.2 GB
NARROW=1 python3 driver.py 16   # 4-byte residues: half the memory, roughly twice the time
```

Useful environment variables:

| variable | effect |
|---|---|
| `NARROW=1` | 30-bit primes and 4-byte residues instead of 61-bit and 8-byte. Halves memory, roughly doubles the number of primes needed. |
| `BIN=./constp` | Use the flat-array solver instead of the layered one. Faster, about four times the memory. |
| `PROGRESS=1` | Report each dynamic programming layer as it finishes. |
| `OMP_NUM_THREADS` | Threads for the layered solver. |

To see what a given `n` will cost before committing to it:

```
$ ./constp2 16 --report
n=16  peak_states=601080390  wide_GB=4.48  narrow_GB=2.24  dense_GB=17.18
```

## Layout

```
constp.cpp           flat array over all 2^(2n-1) subset masks, one prime at a time
constp2.cpp          the same recursion layered by popcount, with colex ranking
driver.py            picks primes, runs a solver, does the CRT, prints exact integers
reference.py         the definition transcribed directly, slow, used as ground truth
verify.py            test suite: solvers against each other, against reference, against data
identity_check.py    symbolic expansion of the operator identity for n = 1, 2
data/const.txt       computed values
docs/algorithm.md    what changed relative to a straightforward memoised implementation
docs/equivalence.md  line-by-line argument that the fast version computes the same thing
cluster/             Slurm job array, one prime per task, and the CRT combiner
tn-parity/           a separate line of work on the mod p reduction, see below
```

## How it works

The starting point is a memoised recursion over the set of steps not yet used. Write
`V = {-(n-1), ..., n-1}` for the step values and `M` for the subset still available:

```
b[{}] = 1
b[M]  = sum over i in M with ns = s(M) + v_i >= 0  of  (-1)^#{j in M : j > i} * g(ns) * b[M \ {i}]

s(M) = -sum_{i in M} v_i          the current height
g(a) = (a+1)(a+2)...(a+n)         the weight
raw(n) = b[V],   const(n) = |raw(n)| / prod_{k<2n} k!
```

Written straightforwardly, in Python or Julia, this runs out of memory around `n = 13`.
The ceiling is not the number of states, it is the cost of each one. Three changes,
none of which touch the mathematics:

**The key is one machine word.** The height `s` and the depth `d` are both determined by
`M`, since the step values sum to zero, so the natural memo key `(avail, s, d)` carries
two redundant fields. Measured at `n = 12`, the keys cost 443 MB against 330 MB for the
values they index, and a tuple or vector key allocates on every branch. As a bitmask the
state is a single `uint64` for every `n` up to 32. Every submask is numerically smaller
than its superset, so filling a flat array in increasing numeric order reaches every
dependency before its dependent, and no hash table is needed at all.

**Two layers, not the whole table.** `b[M]` only ever reads masks of popcount `|M| - 1`,
so two consecutive layers suffice. Addressing inside a layer uses the colex rank of a
`k`-mask with bit positions `p_1 < ... < p_k`, which is `sum_j C(p_j, j)`. Deleting
`p_j` gives `sum_{i<j} C(p_i, i) + sum_{i>j} C(p_i, i-1)`, a prefix sum plus a suffix
sum, so one `O(k)` pass per mask makes all `k` of its transitions `O(1)`. Peak memory
drops from `2^(2n-1)` to `C(2n-1, n) + C(2n-1, n-1)`.

**Residues, not big integers.** The recursion uses only addition, subtraction and
multiplication, and the final division is by a known constant, so the whole thing runs
in `Z/qZ`. That is 8 bytes per state instead of a few hundred, with `__int128` multiplies
and no allocation in the inner loop, and the exact integer is reconstructed by CRT at the
end. Reconstructing `const` rather than `raw` roughly halves the number of primes needed,
since `raw` carries all 1274 bits of `W = prod_{k<2n} k!` as dead weight at `n = 15` and
`W` is invertible modulo any `q > 2n`. The primes are completely independent of one
another, which is where the parallelism lives.

Worth stating plainly: there is no sparsity to exploit here. The weight `(a+1)...(a+n)`
is strictly positive, so nothing cancels, and the number of live states comes out at
`0.52 * 2^(2n-1)` for every `n` tested. The win has to come from making each state cheap,
not from having fewer of them.

`docs/algorithm.md` covers this in more detail, and `docs/equivalence.md` argues step by
step that the result is unchanged.

## Memory and time

Measured on two cores. "Straightforward" is the entry count times the average value size
plus the key words, measured exactly up to `n = 12` and extrapolated above. It ignores
allocator and hash table overhead, so it is an underestimate.

| n | straightforward | `constp` | `constp2` | `constp2` narrow | wall clock |
|---|---|---|---|---|---|
| 13 | ~6 GB | 0.26 GB | 0.09 GB | 0.045 GB | 62 s |
| 14 | ~25 GB | 1.0 GB | 0.30 GB | 0.15 GB | 199 s |
| 15 | ~100 GB | 4.0 GB | 1.16 GB | 0.58 GB | ~12 min |
| 16 | ~450 GB | 16 GB | 4.48 GB | 2.24 GB | ~50 min |
| 17 | | 64 GB | 17.4 GB | 8.7 GB | ~3 h |
| 18 | | 256 GB | 67.6 GB | 33.8 GB | ~12 h |

A 16 GB laptop reaches `n = 16`. A 128 GB cluster node reaches `n = 18`.

## Correctness

The number of primes is sized from a proved upper bound on `|const(n)|` rather than from
a heuristic, so the reconstruction is exact and not merely probable. The bound is
`cluster/primes.py:bound_bits_const`, and `docs/equivalence.md` gives the derivation.

`python3 verify.py` runs four independent comparisons:

- `constp` against `constp2`. The two solvers were written separately and agree residue
  for residue at every `n` and every prime tested, and at 1, 2, 8 and 16 threads.
- Both against `reference.py`, which is the definition transcribed directly, with none of
  the representation changes in it.
- Both against `data/const.txt`, whose first thirteen entries are the OEIS b-file.
- On the cluster path, `cluster/combine.py` re-runs the CRT with the last two primes
  dropped and refuses to print anything if the answer moves.

`identity_check.py` closes the loop at the other end by expanding the operator identity
symbolically for `n = 1` and `n = 2`, so the recursion is checked against the object it
is supposed to be counting and not only against itself.

## Running on a cluster

The primes are independent, so one prime per array task is the entire parallelisation
strategy. `cluster/` holds a Slurm job array, a prime list generator that prints the
array size to use, and a combiner that refuses to guess when tasks are missing. Details
in [`cluster/README.md`](cluster/README.md).

```
bash cluster/build.sh                          # once, on a login node
python3 cluster/primes.py 17 > primes_17.txt   # prints the --array= size to use
mkdir -p logs out/n17
sbatch --array=1-41 cluster/task.sh 17 primes_17.txt out/n17
python3 cluster/combine.py 17 out/n17 primes_17.txt
```

## `tn-parity/`

A separate, self-contained line of work, not needed to compute `const(n)`. It concerns
the reduced sum

```
T_n = sum_{sigma} sign(sigma) prod_{k=1}^{2n-1} floor(E_k(sigma)/n)
```

which controls the `p`-adic valuation of `const(p)` at the prime `p`. Its support turns
out to be a regular language on four states, which makes the states rankable and brings
`n = 23` down from roughly 150 GB to 0.71 GB. Because the resulting value is closed form,
there is also a certificate that runs in `O(n)` memory and a few kilobytes of working
set. See [`tn-parity/README.md`](tn-parity/README.md).

## Citing

If this is useful, please cite the paper:

> K. C. Shah and A. V. Kiselev, *The alternating compositions of weighted differential
> operators yield the weights' Wronskian with which constant?*, arXiv:2605.11137
> [math.CO], 2026.

and, for the sequence itself, [OEIS A392714](https://oeis.org/A392714). Machine readable
metadata is in `CITATION.cff`.

## License

MIT. See [LICENSE](LICENSE).
