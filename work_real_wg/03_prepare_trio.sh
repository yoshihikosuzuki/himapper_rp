#!/usr/bin/env bash
#SBATCH -J hg002_wg_trio
#SBATCH -o 03_prepare_trio.log
#SBATCH -p rl8
#SBATCH -c 128
#SBATCH -n 1
#SBATCH -N 1
set -euo pipefail
source config.sh
module load ${_YAK}
set -x
THREADS=${SLURM_CPUS_PER_TASK:-128}
TRIOBIN="$D_TRIO/$READS_BASE.triobin.tsv"

run_timed "$D_TRIO/triobin.time" \
  yak triobin -t"$THREADS" "$D_TRIO/pat.yak" "$D_TRIO/mat.yak" "$READS" \
  > "$TRIOBIN"

awk 'BEGIN{OFS="\t"}
 NR==FNR {map_hap[$1]=$2; next}
 !($1 in map_hap) {next}
 {h=$5+$6
  if(h==0) s="no_hapmer"
  else if($2!="p" && $2!="m") s="ambiguous"
  else if(($2=="p"?"pat":"mat")!=map_hap[$1]) s="disagree"
  else s=(h<=40?"vhard":(h<=100?"hard":(h<=350?"med":"easy")))
  print $1,s}' "$READS_CONFIDENT_MAPPING_TRUTH" "$TRIOBIN" \
  > "$D_EVAL/hapmer_strata.tsv"
