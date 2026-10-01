#!/usr/bin/env bash
#SBATCH -J hg002_wg_evaluate
#SBATCH -o 08_evaluate.%j.log
#SBATCH -p rl8
#SBATCH -c 4
#SBATCH -n 1
#SBATCH -N 1
set -euo pipefail
source config.sh
module load "${_PYTHON}"
set -x

TOOLS=(
  himapper="$D_HIMAPPER/ref-seed_read-seed/mapping.tsv"
  minimap2="$D_MM2/mapping.tsv"
  wm="$D_WM/mapping.tsv"
  wm2="$D_WM2/mapping.tsv"
)

RESOURCE_ARGS=()
FUSED_HIMAPPER_TIME="$D_HIMAPPER/ref-seed_read-seed/seed-map.time"
if [[ -s "$FUSED_HIMAPPER_TIME" ]]; then
  RESOURCE_ARGS+=(
    --fused-read-seed-map-time "himapper=$FUSED_HIMAPPER_TIME"
  )
fi

SUMMARY_TMP="$D_EVAL/.summary.${SLURM_JOB_ID:-$$}.tmp"
trap 'rm -f "$SUMMARY_TMP"' EXIT
python3 "$SCRIPTS/eval.py" \
  --mapping-truth "$READS_CONFIDENT_MAPPING_TRUTH" \
  --exclude-truth-bed "$EVAL_EXCLUSION_BED" \
  --strata hapmer \
  --strata alignment \
  --strata-file "hapmer=$D_EVAL/hapmer_strata.tsv" \
  --strata-file "alignment=$EVAL_ALIGNMENT_STRATA" \
  --liftover-paf "$D_REF/mat_vs_pat.paf" \
  --liftover-paf "$D_REF/pat_vs_mat.paf" \
  --time-dir "$WR" \
  "${RESOURCE_ARGS[@]}" \
  "${TOOLS[@]}" \
  | tee "$SUMMARY_TMP"
mv "$SUMMARY_TMP" "$D_EVAL/summary.txt"
trap - EXIT
