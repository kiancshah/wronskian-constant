#!/bin/bash
#SBATCH --job-name=constp
#SBATCH --partition=regular
#SBATCH --output=logs/constp_%A_%a.out
#SBATCH --error=logs/constp_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=08:00:00
#SBATCH --mem=24G
#
# One prime per array task. The primes are independent, so this is the parallelism.
# Set --mem from `./constp2 <p> --report`; see SUBMIT.md for the table.
set -e
cd "$SLURM_SUBMIT_DIR"
P=$1; PLIST=$2; OUTDIR=$3
Q=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$PLIST")
[ -n "$Q" ] || { echo "no prime on line $SLURM_ARRAY_TASK_ID of $PLIST"; exit 2; }
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export PROGRESS=1
[ -z "$NARROW" ] && export WIDE=1
./constp2 "$P" "$Q" > "$OUTDIR/q${Q}.txt"
echo "task $SLURM_ARRAY_TASK_ID done: p=$P q=$Q"
