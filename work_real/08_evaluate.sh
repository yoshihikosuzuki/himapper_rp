#!/usr/bin/env bash
#SBATCH -J hg002_evaluate
#SBATCH -o 08_evaluate.log
#SBATCH -p rl8
#SBATCH -c 4
#SBATCH -n 1
#SBATCH -N 1
# Stage 8 — compare the three Himapper seed configurations and both external
# baselines.  Mapping input is the common chr6 read set; scoring is restricted
# to the assigned-locus/resolvable population.
set -euo pipefail
source config.sh
THREADS=${SLURM_CPUS_PER_TASK:-4}

LOCI="$READS_CONFIDENT_LOCUS"
RES="$D_EVAL/resolvability_strata.tsv"
STR="$D_EVAL/hapmer_strata.tsv"
TOOLS=(
  ref-all_read-seed="$D_HIMAPPER/ref-all_read-seed/placement.tsv"
  ref-seed_read-all="$D_HIMAPPER/ref-seed_read-all/placement.tsv"
  ref-seed_read-seed="$D_HIMAPPER/ref-seed_read-seed/placement.tsv"
  mm2="$D_MM2/placement.tsv"
  wm2="$D_WM2/placement.tsv"
)

strata_opt=(--strata "$STR")

# (1) strict-resolvable population
python3 "$SCRIPTS/eval.py" \
  --loci "$LOCI" --resolvable-mode kmer --resolvable "$RES" "${strata_opt[@]}" \
  --liftover-paf "$D_REF/mat_vs_pat.paf" --liftover-paf "$D_REF/pat_vs_mat.paf" \
  --ov 0.5 --time-dir "$WR" "${TOOLS[@]}" \
  | tee "$D_EVAL/stratified.txt"

# (2) variant-containing population
VAR="$D_REF/variants.tsv"
{ echo; echo "########## alignment-resolvable population ##########"; } \
  | tee -a "$D_EVAL/stratified.txt"
python3 "$SCRIPTS/eval.py" \
  --loci "$LOCI" --resolvable-mode alignment --variants "$VAR" "${strata_opt[@]}" \
  --liftover-paf "$D_REF/mat_vs_pat.paf" --liftover-paf "$D_REF/pat_vs_mat.paf" \
  --ov 0.5 --time-dir "$WR" "${TOOLS[@]}" \
  | tee -a "$D_EVAL/stratified.txt"
