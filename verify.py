"""Test suite.

Three independent checks, in increasing order of cost:

  1. Both solvers against the stored values in `data/const.txt`, which are the terms
     of OEIS A392714 plus the values computed since.
  2. `constp` (flat array over masks) against `constp2` (layered, colex-ranked). The
     two were written independently and must agree residue for residue.
  3. Both against `reference.py`, the direct transcription of the definition.

Usage:

    python3 verify.py            # checks 1 and 2 up to n = 11, check 3 up to n = 9
    python3 verify.py 13         # same, up to n = 13 (slow: check 3 dominates)
    python3 verify.py 14 --fast  # skip check 3

Exit status is 0 if everything passed, 1 otherwise.
"""

import os
import subprocess
import sys
import time

import driver

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data", "const.txt")

# reference.py grows fast enough that this is a sensible ceiling for a routine run
REFERENCE_CEILING = 9


def load_known():
    known = {}
    with open(DATA) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            n, value = line.split()
            known[int(n)] = int(value)
    return known


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    top = int(args[0]) if args else 11
    fast = "--fast" in sys.argv

    known = load_known()
    failures = []

    print(f"checking n = 1..{top}\n")
    print(f"{'n':>3}  {'constp2':>9}  {'constp':>9}  {'reference':>10}  {'vs data':>8}")
    print("-" * 48)

    for n in range(1, top + 1):
        row = [f"{n:3d}"]

        t0 = time.time()
        layered, _, _ = driver.run(n)
        row.append(f"{time.time() - t0:8.2f}s")

        env = dict(os.environ, BIN="./constp")
        t0 = time.time()
        flat, _, _ = driver.run(n, env=env)
        row.append(f"{time.time() - t0:8.2f}s")
        if flat != layered:
            failures.append(f"n={n}: constp and constp2 disagree")

        if not fast and n <= REFERENCE_CEILING:
            import reference
            t0 = time.time()
            ref = reference.const(n)
            row.append(f"{time.time() - t0:9.2f}s")
            if ref != layered:
                failures.append(f"n={n}: reference disagrees with constp2")
        else:
            row.append(f"{'skipped':>10}")

        if n in known:
            row.append(f"{'ok' if known[n] == layered else 'MISMATCH':>8}")
            if known[n] != layered:
                failures.append(f"n={n}: disagrees with data/const.txt")
        else:
            row.append(f"{'new':>8}")

        print("  ".join(row), flush=True)

    print()
    if failures:
        for f in failures:
            print("FAIL:", f)
        return 1
    print(f"all checks passed for n = 1..{top}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
