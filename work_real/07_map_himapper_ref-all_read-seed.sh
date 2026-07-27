#!/usr/bin/env bash
#SBATCH -J ref-all_read-seed
#SBATCH -o 07_map_himapper_ref-all_read-seed.log
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
RUN_DIR="$D_HIMAPPER/ref-all_read-seed"
mkdir -p "$RUN_DIR"
cd "$RUN_DIR"

ln -sf "$D_HIMAPPER"/{ref.dam,.ref.bps,.ref.hdr,.ref.idx} .
ln -sf "$D_HIMAPPER"/{reads.dam,.reads.bps,.reads.hdr,.reads.idx} .

NB=$(awk -F'=' '/blocks/{gsub(/ /,"",$2);print $2}' reads.dam)
export CN_UNIQ_W=11.5 CN_TIER_SRC=1 CN_HASH_TIE=1 CN_DUMP=0

run_timed "$RUN_DIR/map.time" \
  "$HIM" -1 -v -T"$THREADS" \
  ref.dam all \
  reads.dam "$D_HIMAPPER/.reads.minimax.kpos"
: > "$RUN_DIR/mapping.tsv"
for k in $(seq 1 "$NB"); do
  LAshow reads ref.dam "reads.$k.ref.las" 2>/dev/null \
    | python3 "$SCRIPTS/lashow_to_placement.py" "$D_HIMAPPER/reads.names" "$D_HIMAPPER/ref.names" \
    >> "$RUN_DIR/mapping.tsv"
done

plot_coverage "$RUN_DIR/coverage.png" \
  --xlo-mb 54 --xhi-mb 68 \
  minimap2="$D_MM2/mapping.tsv" \
  winnowmap2="$D_WM2/mapping.tsv" \
  himapper="$RUN_DIR/mapping.tsv"
