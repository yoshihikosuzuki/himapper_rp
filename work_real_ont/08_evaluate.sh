#!/usr/bin/env bash
#SBATCH -J hg002_ont_evaluate
#SBATCH -o 08_evaluate.log
#SBATCH -p rl8
#SBATCH -c 4
#SBATCH -n 1
#SBATCH -N 1
set -euo pipefail
source config.sh
module load ${_PYTHON}
set -x

TOOLS=(
  himapper="$D_HIMAPPER/ref-seed_read-seed/mapping.tsv"
  minimap2="$D_MM2/mapping.tsv"
)

python3 "$SCRIPTS/eval.py" \
  --mapping-truth "$READS_CONFIDENT_MAPPING_TRUTH" \
  --strata hapmer --strata alignment \
  --liftover-paf "$D_REF/mat_vs_pat.paf" \
  --liftover-paf "$D_REF/pat_vs_mat.paf" \
  --time-dir "$WR" \
  "${TOOLS[@]}" \
  | tee "$D_EVAL/summary.txt"
