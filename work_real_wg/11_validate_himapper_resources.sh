#!/usr/bin/env bash
#SBATCH -J hg002_wg_resource_validate
#SBATCH -o 11_validate_himapper_resources.%j.log
#SBATCH -p rl8
#SBATCH -c 16
#SBATCH -n 1
#SBATCH -N 1
set -euo pipefail
source config.sh

THREADS=${SLURM_CPUS_PER_TASK:-16}
RUN_TAG=${HIMAPPER_RESOURCE_RUN_TAG:?set HIMAPPER_RESOURCE_RUN_TAG}
RUN_ROOT="$WR/experiments/resource_optimization_w100_t16/$RUN_TAG"
PREP_DIR="$RUN_ROOT/prep"
MAP_DIR="$RUN_ROOT/map"
VALID_DIR="$RUN_ROOT/validation"
BIN_DIR="$RUN_ROOT/bin"
JOB_SCRATCH="$WG_TMP_ROOT/resource_optimization_w100_t16/$RUN_TAG"
DATA_DIR="$JOB_SCRATCH/data"

[[ -s "$PREP_DIR/completed.txt" && -s "$MAP_DIR/completed.txt" ]] || {
  echo "Prep or map completion marker is missing" >&2
  exit 1
}
if [[ -e "$VALID_DIR" ]]; then
  echo "Refusing to reuse validation state: $VALID_DIR" >&2
  exit 1
fi
mkdir -p "$VALID_DIR"

set -x
"$TIME_BIN" -v -o "$VALID_DIR/fastk-artifacts-compare.time" \
  bash -c '
    set -euo pipefail
    baseline=$1
    data=$2
    threads=$3
    for leaf in reads.hist reads.ktab reads.prof reads.ref.prof; do
      cmp "$baseline/$leaf" "$data/$leaf"
    done
    cmp "$baseline/ref.reads.prof" "$data/ref.reads.prof"
    for part in $(seq 1 "$threads"); do
      for leaf in \
        ".reads.ktab.$part" \
        ".reads.pidx.$part" \
        ".reads.prof.$part" \
        ".reads.ref.pidx.$part" \
        ".reads.ref.prof.$part" \
        ".ref.reads.pidx.$part" \
        ".ref.reads.prof.$part"; do
        cmp "$baseline/$leaf" "$data/$leaf"
      done
    done
  ' _ "$D_HIMAPPER" "$DATA_DIR" "$THREADS"

"$TIME_BIN" -v -o "$VALID_DIR/ref-profile-decode-compare.time" \
  "$BIN_DIR/compare_profiles" "$D_HIMAPPER/ref.reads" "$DATA_DIR/ref.reads"
"$TIME_BIN" -v -o "$VALID_DIR/reads-profile-decode-compare.time" \
  "$BIN_DIR/compare_profiles" "$D_HIMAPPER/reads" "$DATA_DIR/reads"
"$TIME_BIN" -v -o "$VALID_DIR/reads-ref-profile-decode-compare.time" \
  "$BIN_DIR/compare_profiles" "$D_HIMAPPER/reads.ref" "$DATA_DIR/reads.ref"

{
  printf '%s  %s\n' \
    f66e783a17fbdce78b6b0939c9e20107a24a783b1fb32edb569e05105fe2b1a6 \
    "$DATA_DIR/.ref.minimax.kpos"
  printf '%s  %s\n' \
    62ec55ea876d723ec64266047a4e39612f1fdf5970d4cdfa5a4fdf6411dc77db \
    "$DATA_DIR/.ref.minimax.kpos.idx"
  printf '%s  %s\n' \
    534d9784a1ecaefcf6f376fc1ea3efb220b3dc6c420256c9c90defb168d939c5 \
    "$DATA_DIR/.reads.minimax.kpos"
  printf '%s  %s\n' \
    69da6b7374d7a04b8e5add954dca9253c8cc2fceb15c0e77b200fd74dd0c709f \
    "$DATA_DIR/.reads.minimax.kpos.idx"
  printf '%s  %s\n' \
    80da0a13699e4fc5409624cd1fe77d80be40e3497f26713c6ab1059559678b00 \
    "$MAP_DIR/reads.ref.las"
} > "$VALID_DIR/expected.sha256"

"$TIME_BIN" -v -o "$VALID_DIR/output-sha256.time" \
  sha256sum --check "$VALID_DIR/expected.sha256"
set +x

printf 'FASTK_ARTIFACTS_BYTE_EXACT=PASS\n' > "$VALID_DIR/result.txt"
printf 'PROFILES_FULL_DECODE_EXACT=PASS\n' >> "$VALID_DIR/result.txt"
printf 'KPOS_AND_LAS_SHA256_EXACT=PASS\n' >> "$VALID_DIR/result.txt"
printf 'completed_at=%s\njob_id=%s\n' "$(date --iso-8601=seconds)" \
  "$SLURM_JOB_ID" > "$VALID_DIR/completed.txt"
