#!/usr/bin/env bash
#SBATCH -J hg002_evaluate
#SBATCH -o 08_evaluate.log
#SBATCH -p rl8
#SBATCH -c 4
#SBATCH -n 1
#SBATCH -N 1
set -euo pipefail
source config.sh
set -x

TOOLS=(
#   ref-all_read-seed="$D_HIMAPPER/ref-all_read-seed/mapping.tsv"
  ref-seed_read-seed="$D_HIMAPPER/ref-seed_read-seed/mapping.tsv"
  mm2="$D_MM2/mapping.tsv"
  wm="$D_WM/mapping.tsv"
  wm2="$D_WM2/mapping.tsv"
)

# Score the complete high-confidence mapping-truth population once, then report
# both independent hapmer evidence and alignment self-map difficulty strata.
python3 "$SCRIPTS/eval.py" \
  --mapping-truth "$READS_CONFIDENT_MAPPING_TRUTH" \
  --exclude-truth-bed "$REPO/work_real_wg/eval/hg002v1.1.rDNA_GAP_and_active_HOR.bed" \
  --strata hapmer --strata alignment \
  --liftover-paf "$D_REF/mat_vs_pat.paf" --liftover-paf "$D_REF/pat_vs_mat.paf" \
  "${TOOLS[@]}" \
  | tee "$D_EVAL/summary.txt"
