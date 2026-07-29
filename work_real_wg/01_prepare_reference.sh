#!/usr/bin/env bash
#SBATCH -J hg002_wg_ref
#SBATCH -o 01_prepare_reference.log
#SBATCH -p rl8
#SBATCH -c 16
#SBATCH -n 1
#SBATCH -N 1
set -euo pipefail
source config.sh
module load ${_MINIMAP2}
set -x
THREADS=${SLURM_CPUS_PER_TASK:-16}

# Query MAT against PAT, then query PAT against MAT. eval.py needs the cg tag.
run_timed "$D_REF/mat-vs-pat.time" \
  minimap2 -I8G -c --secondary=no -x asm5 -t"$THREADS" \
  "$D_REF/hg002.pat.fa.gz" "$D_REF/hg002.mat.fa.gz" \
  > "$D_REF/mat_vs_pat.paf"

run_timed "$D_REF/pat-vs-mat.time" \
  minimap2 -I8G -c --secondary=no -x asm5 -t"$THREADS" \
  "$D_REF/hg002.mat.fa.gz" "$D_REF/hg002.pat.fa.gz" \
  > "$D_REF/pat_vs_mat.paf"
