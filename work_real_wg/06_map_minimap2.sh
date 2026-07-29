#!/usr/bin/env bash
#SBATCH -J hg002_wg_minimap2
#SBATCH -o 06_map_minimap2.log
#SBATCH -p rl8
#SBATCH -c 16
#SBATCH -n 1
#SBATCH -N 1
set -euo pipefail
source config.sh
module load ${_MINIMAP2}
set -x
THREADS=${SLURM_CPUS_PER_TASK:-16}

run_timed "$D_MM2/index.time" \
  minimap2 -x map-hifi -k"$MM2_K" -w"$MM2_W" -I8G -t"$THREADS" \
  -d "$D_MM2/dip.mmi" "$REF"

run_timed "$D_MM2/map.time" \
  minimap2 -x map-hifi -k"$MM2_K" -w"$MM2_W" -t"$THREADS" \
  --secondary=no -c "$D_MM2/dip.mmi" "$READS" \
  > "$D_MM2/map.paf"

awk -f "$SCRIPTS/paf_to_placement.awk" "$D_MM2/map.paf" \
  > "$D_MM2/mapping.tsv"
