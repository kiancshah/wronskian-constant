#!/bin/bash
# Run once, on a Slurm login node, from the directory containing constp2.cpp.
#   bash cluster/build.sh
set -e
cd "$(dirname "$0")/.."

# Lmod module names differ between clusters and change over time. Try the usual
# candidates, then fall back to whatever g++ is on PATH.
module purge 2>/dev/null || true
for m in GCC/12.3.0 GCC/11.3.0 GCC foss/2023a foss/2022a; do
  if module load "$m" 2>/dev/null; then echo "loaded module: $m"; break; fi
done

command -v g++ >/dev/null || { echo "no g++ found. run 'module spider GCC' and load one."; exit 1; }
echo "compiler: $(g++ --version | head -1)"

# x86-64-v2 (SSE4.2/POPCNT) is portable across mixed-generation partitions. The hot loop is
# scalar __int128 arithmetic, so -march=native buys nothing and would risk SIGILL if an
# array task lands on an older node than the one you built on.
g++ -O3 -march=x86-64-v2 -std=c++17 -fopenmp -o constp2 constp2.cpp
g++ -O3 -march=x86-64-v2 -std=c++17            -o constp  constp.cpp
echo "built ./constp2 and ./constp"
echo
echo "memory needed per task (this is RAM, not disk):"
for k in 15 16 17 18; do ./constp2 $k --report; done
