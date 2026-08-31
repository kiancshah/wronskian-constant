#!/usr/bin/env python3
"""Emit the prime list for const(n), one per line, sized from a PROVED bound on |const(n)|.

Fixed list => fixed job-array size, and no probabilistic early-termination step anywhere.
    python3 cluster/primes.py 17            > primes_17.txt   # 61-bit primes (8-byte values)
    python3 cluster/primes.py 17 --narrow   > primes_17.txt   # 30-bit primes (4-byte values)
"""
import sys
from math import factorial, prod

def is_prime(x):
    if x < 2: return False
    for p in (2,3,5,7,11,13,17,19,23,29,31,37):
        if x % p == 0: return x == p
    d, r = x-1, 0
    while d % 2 == 0: d //= 2; r += 1
    for a in (2,3,5,7,11,13,17,19,23,29,31,37):
        y = pow(a, d, x)
        if y in (1, x-1): continue
        for _ in range(r-1):
            y = y*y % x
            if y == x-1: break
        else: return False
    return True

def bound_bits_const(p):
    """Proposition 6.2 of constp-computation.tex, in exact integer arithmetic.

        const(p) <= (N-1)!/W * (Ebar^{falling p})^{N-1},   Ebar = V_p/(N-1),
        V_p = sum_k v^m_k,  v^m_k = (k+1)(2p-k)/2  (the paper's master vector).

    Uses only  delta_k >= 0  (so sum_k E_k <= V_p) and concavity of
    E -> log(E^{falling p}), via Jensen. Returns an upper bound on the number of
    binary digits of const(p)."""
    N = 2*p
    M = N - 1
    V = sum((k+1)*(2*p-k)//2 for k in range(1, N))     # sum of the master vector
    assert V > (p-1)*M, "Ebar must exceed p-1 for the falling factorial to be positive"
    # (Ebar^{falling p})^{M} = prod_i (V - i*M)^M / M^(p*M), all integers
    num = factorial(M) * prod(V - i*M for i in range(p))**M
    den = M**(p*M) * prod(factorial(k) for k in range(N))
    bound = num // den + 1                              # ceiling, so still an upper bound
    return max(8, bound.bit_length() + 1)               # +1 for the balanced representative

def primes_for(n, narrow=False):
    bits = 30 if narrow else 61
    need = bound_bits_const(n)
    out, p = [], (1 << bits) + 1
    while sum(q.bit_length() for q in out) < need:
        if is_prime(p): out.append(p)
        p += 2
    return out

if __name__ == "__main__":
    n = int(sys.argv[1]); narrow = "--narrow" in sys.argv
    ps = primes_for(n, narrow)
    sys.stderr.write(f"n={n}: proved bound {bound_bits_const(n)} bits -> {len(ps)} primes "
                     f"({'30' if narrow else '61'}-bit). Submit --array=1-{len(ps)}\n")
    print("\n".join(map(str, ps)))
