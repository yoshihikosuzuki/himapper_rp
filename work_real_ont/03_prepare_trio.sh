#!/usr/bin/env bash
#SBATCH -J hg002_ont_trio
#SBATCH -o 03_prepare_trio.log
#SBATCH -p rl8
#SBATCH -c 128
#SBATCH -n 1
#SBATCH -N 1
set -euo pipefail
source config.sh
module load ${_YAK}
set -x
THREADS=${SLURM_CPUS_PER_TASK:-128}

run_timed "$D_TRIO/triobin.time" \
  yak triobin -t"$THREADS" "$D_TRIO/pat.yak" "$D_TRIO/mat.yak" "$READS" \
  > "$D_TRIO/$READS_BASE.triobin.tsv"
