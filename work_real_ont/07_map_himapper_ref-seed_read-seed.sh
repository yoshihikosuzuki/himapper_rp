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
ln -sf "$D_HIMAPPER"/{ref.dam,.ref.bps,.ref.hdr,.ref.idx} .
ln -sf "$D_HIMAPPER"/{reads.dam,.reads.bps,.reads.hdr,.reads.idx} .
NB=$(awk -F'=' '/blocks/{gsub(/ /,"",$2);print $2}' reads.dam)
REF_LAYOUT="$RUN_DIR/ref.dam.layout"
"$DAZZ_DB_ROOT/DBshow" -n ref.dam > "$REF_LAYOUT"
export CN_UNIQ_W CN_TIER_SRC CN_HASH_TIE CN_DUMP

run_timed "$RUN_DIR/map.time" \
  "$HIM" -k"$KMER" -e"$HIMAPPER_E" -s"$HIMAPPER_TRACE_SPACING" \
  -1 -v -T"$THREADS" -P"$HIMAPPER_TMP" \
  ref.dam "$D_HIMAPPER/.ref.minimax.kpos" \
  reads.dam "$D_HIMAPPER/.reads.minimax.kpos"

LAcat "reads.@1-${NB}.ref.las" > "$RUN_DIR/reads.ref.las"
LAshow reads ref.dam "$RUN_DIR/reads.ref.las" 2>/dev/null \
  | python3 "$SCRIPTS/lashow_to_placement.py" \
      --ref-layout "$REF_LAYOUT" \
      "$D_HIMAPPER/reads.names" "$D_HIMAPPER/ref.names" \
  > "$RUN_DIR/mapping.tsv"
