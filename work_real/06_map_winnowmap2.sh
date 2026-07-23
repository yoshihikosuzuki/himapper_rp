#!/usr/bin/env bash
#SBATCH -J hg002_winnowmap2
#SBATCH -o 06_map_winnowmap2.log
#SBATCH -p rl8
#SBATCH -c 16
#SBATCH -n 1
#SBATCH -N 1
set -euo pipefail
source config.sh
module load ${_WINNOWMAP} ${_SAMTOOLS} ${_MERYL}
set -x
THREADS=${SLURM_CPUS_PER_TASK:-16}
SORT_THREADS=8
BAM_MAP_THREADS=$((THREADS - SORT_THREADS))
(( BAM_MAP_THREADS > 0 )) || {
  echo "SLURM_CPUS_PER_TASK must be greater than SORT_THREADS=$SORT_THREADS" >&2
  exit 1
}

run_timed "$D_WM2/ref-meryl.time" meryl count k=15 output "$D_WM2/chr6.meryl" "$REF"
run_timed "$D_WM2/repetitive-kmers.time" meryl print greater-than distinct=0.9998 "$D_WM2/chr6.meryl" > "$D_WM2/chr6.repetitive.txt"

run_timed "$D_WM2/map.time" \
  winnowmap -W "$D_WM2/chr6.repetitive.txt" -x map-pb -t"$THREADS" --secondary=no -c \
  "$REF" "$READS" > "$D_WM2/map.paf"
awk -f "$SCRIPTS/paf_to_placement.awk" "$D_WM2/map.paf" > "$D_WM2/placement.tsv"

run_timed "$D_WM2/bam.time" \
  bash -o pipefail -c '
    set -x
    winnowmap -W "$1" -x map-pb -t"$2" --secondary=no -a "$3" "$4" \
      | samtools sort -@"$5" -o "$6" -
  ' _ "$D_WM2/chr6.repetitive.txt" "$BAM_MAP_THREADS" "$REF" "$READS" \
    "$SORT_THREADS" "$D_WM2/wm2.sorted.bam"
