#!/usr/bin/env bash
#SBATCH -J hg002_ont_himapper
#SBATCH -o 07_map_himapper_ref-seed_read-seed.log
#SBATCH -p rl8
#SBATCH -c 16
#SBATCH -n 1
#SBATCH -N 1
set -euo pipefail
source config.sh
module load ${_DALIGNER} ${_PYTHON}
set -x
THREADS=${SLURM_CPUS_PER_TASK:-16}
HIM="$REPO/himapper"
RUN_DIR="$D_HIMAPPER/ref-seed_read-seed"
HIMAPPER_TMP="$ONT_TMP_ROOT/himapper"
mkdir -p "$RUN_DIR" "$HIMAPPER_TMP"
cd "$RUN_DIR"

sha256sum "$HIM" "$DAZZ_DB_ROOT/DBshow" \
  "$SCRIPTS/lashow_to_placement.py"
HIM_SHA256=$(sha256sum "$HIM" | awk '{print $1}')
REF_KPOS_SHA256=$(sha256sum "$D_HIMAPPER/.ref.minimax.kpos" \
  | awk '{print $1}')
READS_KPOS_SHA256=$(sha256sum "$D_HIMAPPER/.reads.minimax.kpos" \
  | awk '{print $1}')
ln -sf "$D_HIMAPPER"/{ref.dam,.ref.bps,.ref.hdr,.ref.idx} .
ln -sf "$D_HIMAPPER"/{reads.dam,.reads.bps,.reads.hdr,.reads.idx} .
NB=$(awk -F'=' '/blocks/{gsub(/ /,"",$2);print $2}' reads.dam)
REF_LAYOUT="$RUN_DIR/ref.dam.layout"
"$DAZZ_DB_ROOT/DBshow" -n ref.dam > "$REF_LAYOUT"
export CN_UNIQ_W CN_TIER_SRC CN_HASH_TIE CN_DUMP

{
  printf 'k\t%s\n' "$KMER"
  printf 'seed_w\t%s\n' "$HIMAPPER_WS"
  printf 'e\t%s\n' "$HIMAPPER_E"
  printf 'trace_spacing\t%s\n' "$HIMAPPER_TRACE_SPACING"
  printf 'CN_UNIQ_W\t%s\n' "$CN_UNIQ_W"
  printf 'CN_TIER_SRC\t%s\n' "$CN_TIER_SRC"
  printf 'CN_HASH_TIE\t%s\n' "$CN_HASH_TIE"
  printf 'CN_DUMP\t%s\n' "$CN_DUMP"
  printf 'threads\t%s\n' "$THREADS"
  printf 'mapper_sha256\t%s\n' "$HIM_SHA256"
  printf 'ref_kpos_sha256\t%s\n' "$REF_KPOS_SHA256"
  printf 'reads_kpos_sha256\t%s\n' "$READS_KPOS_SHA256"
} > "$RUN_DIR/parameters.tsv"

run_timed "$RUN_DIR/map.time" \
  "$HIM" -k"$KMER" -e"$HIMAPPER_E" -s"$HIMAPPER_TRACE_SPACING" \
  -1 -v -T"$THREADS" -P"$HIMAPPER_TMP" \
  ref.dam "$D_HIMAPPER/.ref.minimax.kpos" \
  reads.dam "$D_HIMAPPER/.reads.minimax.kpos"

LAcat "reads.@1-${NB}.ref.las" > "$RUN_DIR/reads.ref.las"
# The shared converter is chain-aware by default and emits the maximum
# query-span START/NEXT chain for each read.
LAshow reads ref.dam "$RUN_DIR/reads.ref.las" 2>/dev/null \
  | python3 "$SCRIPTS/lashow_to_placement.py" \
      --ref-layout "$REF_LAYOUT" \
      "$D_HIMAPPER/reads.names" "$D_HIMAPPER/ref.names" \
  > "$RUN_DIR/mapping.tsv"
