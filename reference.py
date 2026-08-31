"""Reference implementation of const(n).

This is the definition, written as directly as it can be written. It is slow and it
uses a lot of memory, and that is the point: it is the thing the fast solvers in
`constp.cpp` and `constp2.cpp` are checked against, so it is kept as literal as
possible rather than optimised.

    const(n) = |raw(n)| / prod_{k=0}^{2n-1} k!

    raw(n)   = sum over orderings sigma of {2, ..., 2n} whose partial sums stay
               admissible, of  sign(sigma) * prod of the falling-factorial weights.

The recursion below is stated over the set `avail` of step values not yet used,
the running height `s`, and the depth `d`. `b(avail, s, d)` is the signed weighted
sum over all completions from that state.

Practical range: n <= 12 or so on a normal machine. Use `driver.py` for anything
larger.
"""

from functools import lru_cache
from math import factorial, prod


def raw(n):
    """The unnormalised alternating sum, before dividing by the Wronskian factor W."""
    N = 2 * n

    @lru_cache(maxsize=None)
    def b(avail, s, d):
        if d == N - 1:
            return 1
        total = 0
        for i, e in enumerate(avail):
            ns = s + e - n - 1                     # new height after taking step e
            if ns < 0:                             # inadmissible, prune
                continue
            sgn = -1 if (len(avail) - 1 - i) & 1 else 1
            total += sgn * prod(range(ns + 1, ns + n + 1)) * b(avail[:i] + avail[i + 1:], ns, d + 1)
        return total

    value = b(tuple(range(2, N + 1)), 0, 0)
    b.cache_clear()
    return value


def W(n):
    """The Wronskian normalisation factor, prod_{k=0}^{2n-1} k!."""
    return prod(factorial(k) for k in range(2 * n))


def const(n):
    """const(n), as an exact Python integer."""
    return abs(raw(n)) // W(n)


if __name__ == "__main__":
    import sys
    top = int(sys.argv[1]) if len(sys.argv) > 1 else 9
    for n in range(1, top + 1):
        print(f"n={n:2d}  const={const(n)}")
