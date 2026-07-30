#!/usr/bin/env bash
#SBATCH -J hg002_wg_resource_prep
#SBATCH -o 09_benchmark_himapper_resources_prep.%j.log
#SBATCH -p rl8
#SBATCH -c 16
#SBATCH -n 1
#SBATCH -N 1
set -euo pipefail
source config.sh

THREADS=${SLURM_CPUS_PER_TASK:-16}
RUN_TAG=${HIMAPPER_RESOURCE_RUN_TAG:?set HIMAPPER_RESOURCE_RUN_TAG}
RUN_PARENT="$WR/experiments/resource_optimization_w100_t16"
RUN_ROOT="$RUN_PARENT/$RUN_TAG"
PREP_DIR="$RUN_ROOT/prep"
BIN_DIR="$RUN_ROOT/bin"
SOURCE_DIR="$RUN_ROOT/source"
BASELINE_DIR="$RUN_ROOT/baseline"
JOB_SCRATCH="$WG_TMP_ROOT/resource_optimization_w100_t16/$RUN_TAG"
DATA_DIR="$JOB_SCRATCH/data"
FASTK_TMP="$JOB_SCRATCH/fastk"

if [[ "$THREADS" != 16 || "$KMER" != 32 || "$HIMAPPER_WS" != 100 ]]; then
  echo "Expected T16, k=32, w=s=100; got T=$THREADS k=$KMER w=$HIMAPPER_WS" >&2
  exit 1
fi
if [[ -e "$RUN_ROOT" || -e "$JOB_SCRATCH" ]]; then
  echo "Refusing to reuse run state: $RUN_ROOT or $JOB_SCRATCH" >&2
  exit 1
fi

check_sha()
{
  local expected=$1
  local path=$2
  local actual

  actual=$(sha256sum "$path" | awk '{print $1}')
  if [[ "$actual" != "$expected" ]]; then
    echo "Unexpected binary SHA-256 for $path: $actual (expected $expected)" >&2
    exit 1
  fi
}

check_sha ee07511eb23b8bdda9ef6083d3d883c3fd54d78e994b1b34c25d015ed553c7fb \
  "$REPO/fastkdual/FastKdual"
check_sha 1d1cfe35a2bda54346caef97d71655dd8fccb6e572062b630972f069a63526f4 \
  "$REPO/fastkdual/compare_profiles"
check_sha 16f62b41a61309c5db9c719a49477d7ccaa4481d0ccdc81a5945f8034b9e3b6a \
  "$REPO/minimaxmer"
check_sha 792c015ed9b0332d9008b039eda91dde80e2e45592ba3e3ad2b975eccaa2bd2f \
  "$REPO/himapper"
check_sha 192f203905d778d486a21995faf2a254c5d9e1459a418f72da7c09c8fdf5b8d8 \
  "$FASTK_ROOT/FastK"

for path in \
  "$D_HIMAPPER/ref.dam" "$D_HIMAPPER/.ref.idx" \
  "$D_HIMAPPER/.ref.bps" "$D_HIMAPPER/.ref.hdr" \
  "$D_HIMAPPER/reads.dam" "$D_HIMAPPER/.reads.idx" \
  "$D_HIMAPPER/.reads.bps" "$D_HIMAPPER/.reads.hdr" \
  "$D_HIMAPPER/ref.ktab" "$D_HIMAPPER/ref.prof" \
  "$D_HIMAPPER/ref.reads.prof" \
  "$D_HIMAPPER/reads.hist" "$D_HIMAPPER/reads.ktab" \
  "$D_HIMAPPER/reads.prof" "$D_HIMAPPER/reads.ref.prof" \
  "$D_HIMAPPER/.ref.minimax.kpos" \
  "$D_HIMAPPER/.ref.minimax.kpos.idx" \
  "$D_HIMAPPER/.reads.minimax.kpos" \
  "$D_HIMAPPER/.reads.minimax.kpos.idx" \
  "$D_HIMAPPER/ref-seed_read-seed/reads.ref.las"; do
  [[ -s "$path" ]] || { echo "Missing baseline input: $path" >&2; exit 1; }
done
for part in $(seq 1 "$THREADS"); do
  for path in \
    "$D_HIMAPPER/.ref.ktab.$part" \
    "$D_HIMAPPER/.ref.pidx.$part" \
    "$D_HIMAPPER/.ref.prof.$part" \
    "$D_HIMAPPER/.ref.reads.pidx.$part" \
    "$D_HIMAPPER/.ref.reads.prof.$part" \
    "$D_HIMAPPER/.reads.ktab.$part" \
    "$D_HIMAPPER/.reads.pidx.$part" \
    "$D_HIMAPPER/.reads.prof.$part" \
    "$D_HIMAPPER/.reads.ref.pidx.$part" \
    "$D_HIMAPPER/.reads.ref.prof.$part"; do
    [[ -s "$path" ]] || { echo "Missing baseline shard: $path" >&2; exit 1; }
  done
done

mkdir -p "$PREP_DIR" "$BIN_DIR" "$SOURCE_DIR/fastkdual" \
  "$SOURCE_DIR/root" "$SOURCE_DIR/jobs" "$BASELINE_DIR" \
  "$DATA_DIR" "$FASTK_TMP"

cp -p "$REPO/fastkdual/FastKdual" "$REPO/fastkdual/compare_profiles" \
  "$REPO/minimaxmer" "$REPO/himapper" "$FASTK_ROOT/FastK" "$BIN_DIR/"
cp -p "$REPO/fastkdual/FastKdual.c" "$REPO/fastkdual/FastKdual.h" \
  "$REPO/fastkdual/count_dual.c" "$REPO/fastkdual/compare_profiles.c" \
  "$REPO/fastkdual/Makefile" "$SOURCE_DIR/fastkdual/"
cp -p "$REPO/minimaxmer.c" "$REPO/map.c" "$REPO/map.h" \
  "$REPO/himapper.c" "$REPO/Makefile" "$SOURCE_DIR/root/"
cp -p "$WR/09_benchmark_himapper_resources_prep.sh" \
  "$WR/10_benchmark_himapper_resources_map.sh" \
  "$WR/11_validate_himapper_resources.sh" \
  "$WR/config.sh" "$SOURCE_DIR/jobs/"
cp -p "$D_HIMAPPER/reads-profile.time" \
  "$D_HIMAPPER/reads-vs-ref-profile.time" \
  "$D_HIMAPPER/ref-vs-reads-profile.time" \
  "$D_HIMAPPER/ref-seeds.time" "$D_HIMAPPER/read-seeds.time" \
  "$D_HIMAPPER/ref-seed_read-seed/map.time" "$BASELINE_DIR/"

sha256sum "$BIN_DIR/FastKdual" "$BIN_DIR/compare_profiles" \
  "$BIN_DIR/minimaxmer" "$BIN_DIR/himapper" "$BIN_DIR/FastK" \
  > "$RUN_ROOT/binaries.sha256"
git -C "$REPO" rev-parse HEAD > "$RUN_ROOT/git-head.txt"
git -C "$REPO" status --short --untracked-files=all > "$RUN_ROOT/git-status.txt"
git -C "$REPO" diff -- . ':(exclude)work_real_wg/experiments' \
  > "$RUN_ROOT/source.diff"
module list > "$RUN_ROOT/modules.txt" 2>&1 || true
{
  printf 'run_tag=%s\n' "$RUN_TAG"
  printf 'prep_job_id=%s\n' "${SLURM_JOB_ID:?}"
  printf 'threads=%s\nkmer=%s\nwindow=%s\nspacing=%s\n' \
    "$THREADS" "$KMER" "$HIMAPPER_WS" "$HIMAPPER_WS"
  printf 'reference_dam=%s\nreads_dam=%s\n' \
    "$D_HIMAPPER/ref.dam" "$D_HIMAPPER/reads.dam"
  printf 'baseline_dir=%s\nscratch=%s\n' "$D_HIMAPPER" "$JOB_SCRATCH"
} > "$RUN_ROOT/run.config"

cd "$PREP_DIR"
set -x
run_timed "$PREP_DIR/reads-dual-profile.time" \
  "$BIN_DIR/FastKdual" -k"$KMER" -v -T"$THREADS" -t1 \
  -p:"$D_HIMAPPER/ref" -N"$DATA_DIR/reads" -R"$DATA_DIR/reads.ref" \
  -P"$FASTK_TMP" "$D_HIMAPPER/reads.dam"

run_timed "$PREP_DIR/ref-vs-reads-profile.time" \
  "$BIN_DIR/FastK" -k"$KMER" -v -T"$THREADS" -t1 \
  -p:"$DATA_DIR/reads" -N"$DATA_DIR/ref.reads" \
  -P"$FASTK_TMP" "$D_HIMAPPER/ref.dam"

run_timed "$PREP_DIR/ref-seeds.time" \
  "$BIN_DIR/minimaxmer" -T"$THREADS" -k"$KMER" \
  -w"$HIMAPPER_WS" -s"$HIMAPPER_WS" \
  "$D_HIMAPPER/ref.dam" "$D_HIMAPPER/ref.prof" \
  "$DATA_DIR/ref.reads.prof" "$DATA_DIR/.ref.minimax.kpos"

run_timed "$PREP_DIR/read-seeds.time" \
  "$BIN_DIR/minimaxmer" -T"$THREADS" -k"$KMER" \
  -w"$HIMAPPER_WS" -s"$HIMAPPER_WS" \
  "$D_HIMAPPER/reads.dam" "$DATA_DIR/reads.ref.prof" \
  "$DATA_DIR/reads.prof" "$DATA_DIR/.reads.minimax.kpos"
set +x

for path in \
  "$DATA_DIR/reads.hist" "$DATA_DIR/reads.ktab" \
  "$DATA_DIR/reads.prof" "$DATA_DIR/reads.ref.prof" \
  "$DATA_DIR/ref.reads.prof" \
  "$DATA_DIR/.ref.minimax.kpos" "$DATA_DIR/.ref.minimax.kpos.idx" \
  "$DATA_DIR/.reads.minimax.kpos" "$DATA_DIR/.reads.minimax.kpos.idx"; do
  [[ -s "$path" ]] || { echo "Missing new prep output: $path" >&2; exit 1; }
done

printf 'completed_at=%s\njob_id=%s\n' "$(date --iso-8601=seconds)" \
  "$SLURM_JOB_ID" > "$PREP_DIR/completed.txt"
