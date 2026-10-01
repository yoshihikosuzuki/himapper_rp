#!/usr/bin/env bash
#SBATCH -J hg002_wg_prepare_eval
#SBATCH -o 04_prepare_eval.%j.log
#SBATCH -p rl8
#SBATCH -c 32
#SBATCH -n 1
#SBATCH -N 1
set -euo pipefail
source config.sh
module load "${_MINIMAP2}" "${_PYTHON}"
set -x

THREADS=${SLURM_CPUS_PER_TASK:-32}
ALIGNMENT_TMP_ROOT="$WG_TMP_ROOT/alignment-strata"
mkdir -p "$ALIGNMENT_TMP_ROOT"

# Rebuild the versioned 0-based, half-open evaluation exclusion mask.
python3 "$REPO/work_real/scripts/build_eval_exclusion_mask.py" \
  --rdna-gap-bed "$EVAL_RDNA_GAP_BED" \
  --censat-bed "$EVAL_CENSAT_BED" \
  --output-bed "$EVAL_EXCLUSION_BED" \
  --metadata "$EVAL_EXCLUSION_METADATA" \
  --reference-label HG002v1.2 \
  --expected-intervals 84 \
  --expected-bp 164865724

# Preserve the chr6 definition: exact truth reference spans self-mapped with
# minimap2 default k/w, secondary candidates retained, near=2 and med=10.
# Generate all Q60 truth reads once; eval.py applies the mask to this stratum.
run_timed "$D_EVAL/alignment-strata.time" \
  python3 "$REPO/work_real/scripts/alignment_strata.py" \
    --threads "$THREADS" \
    --near 2 \
    --med 10 \
    --tmp-dir "$ALIGNMENT_TMP_ROOT" \
    --index-batch-size 8G \
    --minimap2-log "$D_EVAL/alignment-strata.minimap2.log" \
    "$REF" \
    "$READS_CONFIDENT_MAPPING_TRUTH" \
    "$EVAL_ALIGNMENT_STRATA"
