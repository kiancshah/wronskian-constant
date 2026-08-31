"""Single-machine driver: picks the primes, runs the solver, does the CRT.

    python3 driver.py 14                # print const(1) through const(14)
    NARROW=1 python3 driver.py 16       # 4-byte residues: half the memory, twice the primes
    BIN=./constp python3 driver.py 13   # use the flat-array solver instead of the layered one

The solver binaries know nothing about big integers. They are handed one prime at a
time and return raw(n) mod q. Everything below turns those residues into the exact
integer const(n).

The number of primes comes from `cluster/primes.py:bound_bits_const`, a proved upper
bound on the size of const(n), so there is no probabilistic early-termination step
anywhere in this file. Two spare primes are added on top of the bound.
"""

import os
import subprocess
import sys
from math import factorial, prod

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "cluster"))

from primes import bound_bits_const   # single source of truth for the bound


def is_prime(x):
    """Deterministic Miller-Rabin for the range used here."""
    if x < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if x % p == 0:
            return x == p
    d, r = x - 1, 0
    while d % 2 == 0:
        d //= 2
        r += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        y = pow(a, d, x)
        if y in (1, x - 1):
            continue
        for _ in range(r - 1):
            y = y * y % x
            if y == x - 1:
                break
        else:
            return False
    return True


def primes_above(x, cnt):
    """The first `cnt` primes greater than x."""
    out, p = [], x | 1
    while len(out) < cnt:
        if is_prime(p):
            out.append(p)
        p += 2
    return out


def run(n, verbose=True, env=None):
    """Compute const(n) exactly.

    Returns (value, number_of_primes_used, solver_stderr).
    """
    N = 2 * n
    W = prod(factorial(k) for k in range(N))
    need = bound_bits_const(n)

    env = dict(env if env is not None else os.environ)
    narrow = env.get("NARROW") == "1"        # 4-byte values: half the memory, ~2x the primes
    bits = 30 if narrow else 61
    ps = primes_above((1 << bits) + 1, need // bits + 2)
    if not narrow:
        env["WIDE"] = "1"                    # 8-byte values

    binary = env.get("BIN", "./constp2")
    if not os.path.isabs(binary):
        binary = os.path.join(HERE, os.path.basename(binary))
    if not os.path.exists(binary):
        sys.exit(f"{binary} not found. Run `make` first.")

    r = subprocess.run([binary, str(n)] + [str(p) for p in ps],
                       capture_output=True, text=True, env=env)
    if r.returncode:
        print(r.stderr)
        sys.exit(1)

    res = {}
    for line in r.stdout.split("\n"):
        if not line.strip():
            continue
        q, v = line.split()
        q, v = int(q), int(v)
        res[q] = v * pow(W % q, -1, q) % q    # const = raw / W  (mod q)

    M, x = 1, 0
    for q, v in res.items():                  # incremental CRT
        M2 = M * q
        x = (x * q * pow(q, -1, M) + v * M * pow(M, -1, q)) % M2 if M > 1 else v
        M = M2
    if x > M // 2:                            # balanced representative
        x -= M
    return abs(x), len(ps), r.stderr


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    for n in range(1, int(sys.argv[1]) + 1):
        v, k, err = run(n)
        print(f"n={n:2d}  primes={k:2d}  const={v}", flush=True)
