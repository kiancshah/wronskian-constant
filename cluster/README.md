# Running on a Slurm cluster

The primes used in the modular reconstruction are completely independent of one another,
so one prime per array task is the whole parallelisation strategy. There is no
communication and no reduction step until the very end, when the residues are combined
on a login node in a second or two.

These scripts were written for Habrók at the University of Groningen. Nothing in them is
specific to that cluster beyond the module names in `build.sh`, which fall back to
whatever `g++` is on the path.

## 1. Build, once, on a login node

```
bash cluster/build.sh
```

It loads a GCC module, compiles both solvers, and prints the memory each `n` will need.

Two choices worth knowing about:

- Built with `-march=x86-64-v2`, not `-march=native`. The hot loop is scalar `__int128`
  arithmetic, so native instructions buy nothing, and an array task landing on an older
  node than the one you compiled on would take a SIGILL.
- Allocation failure is caught and reported with the layer number and the exact number of
  GB needed, so an undersized `--mem` gives you a message instead of an OOM kill.

## 2. Generate the prime list

```
python3 cluster/primes.py 17 > primes_17.txt
```

The list length comes from a proved bound on `|const(n)|`, so it is fixed in advance,
which is what makes a fixed job array size correct. It prints the array size to use on
stderr, for example `Submit --array=1-41`.

For 4-byte residues, generate the list to match:

```
python3 cluster/primes.py 17 --narrow > primes_17.txt
```

Mixing the two is caught rather than silently truncated: the solver refuses a prime above
32 bits when values are 4 bytes wide.

## 3. Submit

```
mkdir -p logs out/n17
sbatch --array=1-41 cluster/task.sh 17 primes_17.txt out/n17
```

Each task writes one line, `<prime> <raw(n) mod prime>`, into `out/n17/`. Set `--mem`
from the table below, or from `./constp2 <n> --report`.

| n | primes needed | `--mem` wide | `--mem` narrow | wall clock per task, 8 cores |
|---|---|---|---|---|
| 16 | 36 | 8G | 6G | ~5 min |
| 17 | 41 | 24G | 14G | ~15 min |
| 18 | 47 | 80G | 44G | ~50 min |

`NARROW=1 sbatch ...` halves the memory at the cost of roughly twice as many primes.

## 4. Combine

```
python3 cluster/combine.py 17 out/n17 primes_17.txt
```

`combine.py` does not guess. If any task failed it names the exact `--array=` list to
rerun and exits. If everything is present it does the CRT, then re-runs the same
reconstruction with the last two primes dropped and refuses to print an answer if the two
disagree, which would mean either a bad residue or an undersized bound.
