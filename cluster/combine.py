#!/usr/bin/env python3
"""CRT the per-prime residues into const(n).

    python3 cluster/combine.py 17 out/n17 primes_17.txt
Each file in out/n17 is one line: "<prime> <raw(n) mod prime>", as written by task.sh.
"""
import sys, os
from math import factorial, prod

def crt(pairs):
    M, x = 1, 0
    for q, v in pairs:
        M2 = M*q
        x = (x*q*pow(q,-1,M) + v*M*pow(M,-1,q)) % M2 if M > 1 else v % q
        M = M2
    return (x - M if x > M//2 else x), M

n = int(sys.argv[1]); outdir = sys.argv[2]; plist = sys.argv[3]
want = [int(l) for l in open(plist) if l.strip()]
W = prod(factorial(k) for k in range(2*n))

got = {}
for fn in os.listdir(outdir):
    txt = open(os.path.join(outdir, fn)).read().split()
    if len(txt) != 2: continue
    q, v = int(txt[0]), int(txt[1])
    got[q] = v * pow(W % q, -1, q) % q          # const = raw * W^-1  (mod q)

missing = [q for q in want if q not in got]
if missing:
    idx = [want.index(q)+1 for q in missing]
    sys.exit(f"missing {len(missing)} residue(s); rerun --array={','.join(map(str,idx))}")

pairs = [(q, got[q]) for q in want]
x, M = crt(pairs)
x2, _ = crt(pairs[:-2])                          # free cross-check on top of the proved bound
if x != x2:
    sys.exit("FAIL: reconstruction changed when the last two primes were dropped. "
             "A residue is wrong, or the bound was under-sized.")
v = abs(x)
print(f"n={n}  primes={len(pairs)}  digits={len(str(v))}")
print(f"const({n}) = {v}")
