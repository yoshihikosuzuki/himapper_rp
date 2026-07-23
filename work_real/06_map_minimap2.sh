#!/usr/bin/env bash
#SBATCH -J hg002_minimap2
#SBATCH -o 06_map_minimap2.log
#SBATCH -p rl8
#SBATCH -c 16
#SBATCH -n 1
#SBATCH -N 1
set -euo pipefail
source config.sh
module load ${_MINIMAP2} ${_SAMTOOLS}
set -x
THREADS=${SLURM_CPUS_PER_TASK:-16}
SORT_THREADS=8
BAM_MAP_THREADS=$((THREADS - SORT_THREADS))
(( BAM_MAP_THREADS > 0 )) || {
  echo "SLURM_CPUS_PER_TASK must be greater than SORT_THREADS=$SORT_THREADS" >&2
  exit 1
}

run_timed "$D_MM2/index.time" minimap2 -x map-hifi -t"$THREADS" -d "$D_MM2/chr6.mmi" "$REF"

run_timed "$D_MM2/map.time" \
  minimap2 -x map-hifi -t"$THREADS" --secondary=no -c "$D_MM2/chr6.mmi" "$READS" \
  > "$D_MM2/map.paf"
awk -f "$SCRIPTS/paf_to_placement.awk" "$D_MM2/map.paf" > "$D_MM2/placement.tsv"

run_timed "$D_MM2/bam.time" \
  bash -o pipefail -c '
    set -x
    minimap2 -x map-hifi -t"$1" --secondary=no -a "$2" "$3" \
      | samtools sort -@"$4" -o "$5" -
  ' _ "$BAM_MAP_THREADS" "$D_MM2/chr6.mmi" "$READS" \
    "$SORT_THREADS" "$D_MM2/mm2.sorted.bam"
