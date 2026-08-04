#!/usr/bin/env bash
#SBATCH -J ref-seed_read-all
#SBATCH -o 07_map_himapper_ref-seed_read-all.log
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
RUN_DIR="$D_HIMAPPER/ref-seed_read-all"
mkdir -p "$RUN_DIR"
cd "$RUN_DIR"

ln -sf "$D_HIMAPPER"/{ref.dam,.ref.bps,.ref.hdr,.ref.idx} .
ln -sf "$D_HIMAPPER"/{reads.dam,.reads.bps,.reads.hdr,.reads.idx} .

NB=$(awk -F'=' '/blocks/{gsub(/ /,"",$2);print $2}' reads.dam)
REF_LAYOUT="$RUN_DIR/ref.dam.layout"
"$DAZZ_DB_ROOT/DBshow" -n ref.dam > "$REF_LAYOUT"
export CN_TIER_SRC=0 CN_HASH_TIE=1 CN_DUMP=0

run_cfg() {
  local name=$1 weight=$2
  export CN_UNIQ_W="$weight"
  run_timed "$RUN_DIR/map.$name.time" \
    "$HIM" -1 -T"$THREADS" -e"$HIMAPPER_E" \
    ref.dam "$D_HIMAPPER/.ref.minimax.kpos" \
    reads.dam all
  LAcat "reads.@1-${NB}.ref.las" > "$RUN_DIR/reads.ref.$name.las"
  LAshow reads ref.dam "$RUN_DIR/reads.ref.$name.las" 2>/dev/null \
    | python3 "$SCRIPTS/lashow_to_placement.py" \
        --ref-layout "$REF_LAYOUT" \
        "$D_HIMAPPER/reads.names" "$D_HIMAPPER/ref.names" \
    > "$RUN_DIR/mapping.$name.tsv"
  for k in $(seq 1 "$NB"); do
    mv -f "reads.$k.ref.las" "reads.$k.ref.$name.las"
  done
}

run_cfg w8 "$CN_UNIQ_W"

# Stable primary-comparison aliases point to the configured W=8 candidate.
ln -sfn mapping.w8.tsv "$RUN_DIR/mapping.tsv"
ln -sfn map.w8.time "$RUN_DIR/map.time"

plot_coverage "$RUN_DIR/coverage.png" \
  --xlo-mb 54 --xhi-mb 68 \
  minimap2="$D_MM2/mapping.tsv" \
  winnowmap2="$D_WM2/mapping.tsv" \
  w8="$RUN_DIR/mapping.w8.tsv"
