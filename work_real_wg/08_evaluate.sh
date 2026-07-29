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
)

# Rebuild and validate the derived mask so the summary always records the
# current source annotations. Alignment strata cover the full Q60 truth set
# and are filtered to the same masked population inside eval.py.
python3 "$REPO/work_real/scripts/build_eval_exclusion_mask.py" \
  --rdna-gap-bed "$EVAL_RDNA_GAP_BED" \
  --censat-bed "$EVAL_CENSAT_BED" \
  --output-bed "$EVAL_EXCLUSION_BED" \
  --metadata "$EVAL_EXCLUSION_METADATA" \
  --expected-intervals 84 \
  --expected-bp 164859611

for required in \
  "$D_EVAL/hapmer_strata.tsv" \
  "$EVAL_ALIGNMENT_STRATA"; do
  if [[ ! -s "$required" ]]; then
    echo "missing required evaluation stratum: $required" >&2
    exit 2
  fi
done

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
  "${TOOLS[@]}" \
  | tee "$SUMMARY_TMP"
mv "$SUMMARY_TMP" "$D_EVAL/summary.txt"
trap - EXIT
