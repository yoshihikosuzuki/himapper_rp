#!/usr/bin/env bash
#SBATCH -J hg002_wg_resource_map
#SBATCH -o 10_benchmark_himapper_resources_map.%j.log
#SBATCH -p rl8
#SBATCH -c 16
#SBATCH -n 1
#SBATCH -N 1
set -euo pipefail
source config.sh
module load "${_DALIGNER}"

THREADS=${SLURM_CPUS_PER_TASK:-16}
RUN_TAG=${HIMAPPER_RESOURCE_RUN_TAG:?set HIMAPPER_RESOURCE_RUN_TAG}
RUN_ROOT="$WR/experiments/resource_optimization_w100_t16/$RUN_TAG"
PREP_DIR="$RUN_ROOT/prep"
MAP_DIR="$RUN_ROOT/map"
BIN_DIR="$RUN_ROOT/bin"
JOB_SCRATCH="$WG_TMP_ROOT/resource_optimization_w100_t16/$RUN_TAG"
DATA_DIR="$JOB_SCRATCH/data"
HIMAPPER_TMP="$JOB_SCRATCH/himapper"

if [[ "$THREADS" != 16 || "$KMER" != 32 || "$HIMAPPER_WS" != 100 ]]; then
  echo "Expected T16, k=32, w=s=100; got T=$THREADS k=$KMER w=$HIMAPPER_WS" >&2
  exit 1
fi
[[ -s "$PREP_DIR/completed.txt" ]] || {
  echo "Missing completed prep marker: $PREP_DIR/completed.txt" >&2
  exit 1
}
if [[ -e "$MAP_DIR" || -e "$HIMAPPER_TMP" ]]; then
  echo "Refusing to reuse map state: $MAP_DIR or $HIMAPPER_TMP" >&2
  exit 1
fi
for path in \
  "$BIN_DIR/himapper" \
  "$DATA_DIR/.ref.minimax.kpos" "$DATA_DIR/.ref.minimax.kpos.idx" \
  "$DATA_DIR/.reads.minimax.kpos" "$DATA_DIR/.reads.minimax.kpos.idx"; do
  [[ -s "$path" ]] || { echo "Missing map input: $path" >&2; exit 1; }
done

mkdir -p "$MAP_DIR" "$HIMAPPER_TMP"
cd "$MAP_DIR"
ln -s "$D_HIMAPPER"/{ref.dam,.ref.bps,.ref.hdr,.ref.idx} .
ln -s "$D_HIMAPPER"/{reads.dam,.reads.bps,.reads.hdr,.reads.idx} .
NB=$(awk -F'=' '/blocks/{gsub(/ /,"",$2);print $2}' reads.dam)
if [[ "$NB" != 735 ]]; then
  echo "Expected 735 reads blocks, found $NB" >&2
  exit 1
fi

export CN_UNIQ_W CN_TIER_SRC CN_HASH_TIE CN_DUMP
{
  printf 'map_job_id=%s\nthreads=%s\n' "${SLURM_JOB_ID:?}" "$THREADS"
  printf 'CN_UNIQ_W=%s\nCN_TIER_SRC=%s\nCN_HASH_TIE=%s\nCN_DUMP=%s\n' \
    "$CN_UNIQ_W" "$CN_TIER_SRC" "$CN_HASH_TIE" "$CN_DUMP"
  stat -c '%n %s bytes %y' \
    "$DATA_DIR/.ref.minimax.kpos" "$DATA_DIR/.ref.minimax.kpos.idx" \
    "$DATA_DIR/.reads.minimax.kpos" "$DATA_DIR/.reads.minimax.kpos.idx"
} > "$MAP_DIR/inputs.txt"

set -x
run_timed "$MAP_DIR/map.time" \
  "$BIN_DIR/himapper" -k"$KMER" -1 -v -T"$THREADS" -P"$HIMAPPER_TMP" \
  ref.dam "$DATA_DIR/.ref.minimax.kpos" \
  reads.dam "$DATA_DIR/.reads.minimax.kpos"

"$TIME_BIN" -v -o "$MAP_DIR/las-cat.time" \
  LAcat "reads.@1-${NB}.ref.las" > "$MAP_DIR/reads.ref.las"
set +x

[[ -s "$MAP_DIR/reads.ref.las" ]] || {
  echo "Missing consolidated LAS: $MAP_DIR/reads.ref.las" >&2
  exit 1
}
printf 'completed_at=%s\njob_id=%s\n' "$(date --iso-8601=seconds)" \
  "$SLURM_JOB_ID" > "$MAP_DIR/completed.txt"
