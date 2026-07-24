#!/usr/bin/env bash
#SBATCH -J ref-seed_read-all
#SBATCH -o 07_map_himapper_ref-seed_read-all.log
#SBATCH -p rl8
#SBATCH -c 16
#SBATCH -n 1
#SBATCH -N 1
# Himapper reference-seed/read-all W sweep.  W=0 is the no-copy-bonus baseline;
# copy-aware anchoring is always enabled.  W=11.5 is the initial comparison candidate.
set -euo pipefail
source config.sh
module load ${_DALIGNER} ${_PYTHON}
set -x

THREADS=${SLURM_CPUS_PER_TASK:-16}
HIM="$REPO/himapper"
HPC="$REPO/HPC.himapper"
RUN_DIR="$D_HIMAPPER/ref-seed_read-all"
mkdir -p "$RUN_DIR"
cd "$RUN_DIR"

ln -sf "$D_HIMAPPER"/{ref.dam,.ref.bps,.ref.hdr,.ref.idx} .
ln -sf "$D_HIMAPPER"/{reads.dam,.reads.bps,.reads.hdr,.reads.idx} .

NB=$(awk -F'=' '/blocks/{gsub(/ /,"",$2);print $2}' reads.dam)
export CN_TIER_SRC=0 CN_HASH_TIE=1 CN_DUMP=0

COMMANDS="$RUN_DIR/command.sh"
{
  echo '#!/usr/bin/env bash'
  echo 'set -euo pipefail'
  "$HPC" -1 -T"$THREADS" -K"$D_HIMAPPER/.ref.#.minimax.kpos" ref.dam reads 2>/dev/null \
    | grep '^himapper' | sed "s|^himapper|$HIM|"
} > "$COMMANDS"
n_commands=$(tail -n +3 "$COMMANDS" | wc -l)
[[ "$n_commands" -eq "$NB" ]] || {
  echo "generated $n_commands commands for $NB blocks" >&2
  exit 1
}

run_cfg() {
  local name=$1 weight=$2
  export CN_UNIQ_W="$weight"
  run_timed "$RUN_DIR/map.$name.time" bash "$COMMANDS"
  : > "$RUN_DIR/mapping.$name.tsv"
  for k in $(seq 1 "$NB"); do
    LAshow reads ref.dam "reads.$k.ref.las" 2>/dev/null \
      | python3 "$SCRIPTS/lashow_to_placement.py" "$D_HIMAPPER/reads.names" "$D_HIMAPPER/ref.names" \
      >> "$RUN_DIR/mapping.$name.tsv"
    mv -f "reads.$k.ref.las" "reads.$k.ref.$name.las"
  done
}

run_cfg w0    0
run_cfg w1p15 1.15
run_cfg w2p3  2.3
run_cfg w4p6  4.6
run_cfg w11p5 11.5
run_cfg winf  1000000

# Stable primary-comparison aliases point to the explicit W=11.5 candidate.
ln -sfn mapping.w11p5.tsv "$RUN_DIR/mapping.tsv"
ln -sfn map.w11p5.time "$RUN_DIR/map.time"

plot_coverage "$RUN_DIR/coverage.png" \
  --xlo-mb 54 --xhi-mb 68 \
  minimap2="$D_MM2/mapping.tsv" \
  winnowmap2="$D_WM2/mapping.tsv" \
  w0="$RUN_DIR/mapping.w0.tsv" \
  w1p15="$RUN_DIR/mapping.w1p15.tsv" \
  w2p3="$RUN_DIR/mapping.w2p3.tsv" \
  w4p6="$RUN_DIR/mapping.w4p6.tsv" \
  w11p5="$RUN_DIR/mapping.w11p5.tsv" \
  winf="$RUN_DIR/mapping.winf.tsv"
