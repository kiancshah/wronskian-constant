"""Check the operator identity itself, symbolically, for small n.

`driver.py` computes const(n) from a combinatorial recursion. This script does
something different and slower: it expands both sides of the identity that defines
const(n) as actual differential operators and reads the constant off. It is a check
that the recursion is solving the problem it claims to solve, and it doubles as a
worked example of what the number means.

    sum_{sigma in S_{2n}} sign(sigma) w_{sigma(1)}D^n o ... o w_{sigma(2n)}D^n
        = const(n) * Wronskian(w_1, ..., w_{2n}) * D^n

A differential operator is represented as a dict {order: coefficient}, and
composition uses the Leibniz rule

    (a D^i) o (b D^j) = sum_k C(i,k) a b^(k) D^(i+j-k).

Cost is (2n)! compositions of operators of order up to 2n^2, so n = 1 and n = 2 are
comfortable, n = 3 takes a while, and n >= 4 is not worth attempting this way.

    python3 identity_check.py 2

Requires sympy.
"""

import sys
from itertools import permutations
from math import comb

import sympy


def compose(A, B, x):
    """Compose two operators given as {order: coefficient} dicts."""
    out = {}
    for i, a in A.items():
        for j, b in B.items():
            for k in range(i + 1):
                order = i + j - k
                term = comb(i, k) * a * sympy.diff(b, x, k)
                out[order] = out.get(order, 0) + term
    return out


def wronskian(ws, x):
    return sympy.Matrix(len(ws), len(ws),
                        lambda r, c: sympy.diff(ws[c], x, r)).det()


def check(n):
    x = sympy.Symbol("x")
    N = 2 * n
    ws = [sympy.Function(f"w{i + 1}")(x) for i in range(N)]

    total = {}
    for sigma in permutations(range(N)):
        # sign of the permutation
        sign, seen = 1, list(sigma)
        for a in range(N):
            for b in range(a + 1, N):
                if seen[a] > seen[b]:
                    sign = -sign

        term = {n: ws[sigma[0]]}
        for idx in sigma[1:]:
            term = compose(term, {n: ws[idx]}, x)

        for order, coeff in term.items():
            total[order] = total.get(order, 0) + sign * coeff

    # every order other than n must cancel, and order n must be a multiple of the Wronskian
    residue = {o: sympy.simplify(c) for o, c in total.items() if o != n}
    leftovers = {o: c for o, c in residue.items() if c != 0}

    wr = wronskian(ws, x)
    ratio = sympy.simplify(sympy.cancel(sympy.together(total.get(n, 0)) / wr))

    print(f"n = {n}")
    print(f"  orders present besides D^{n}: {sorted(leftovers) if leftovers else 'none'}")
    print(f"  coefficient of D^{n}, divided by the Wronskian: {ratio}")
    return ratio, leftovers


if __name__ == "__main__":
    top = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    ok = True
    for n in range(1, top + 1):
        ratio, leftovers = check(n)
        if leftovers or not ratio.is_Integer:
            ok = False
    print("\nconsistent" if ok else "\nsomething is off")
