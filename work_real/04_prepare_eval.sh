#!/usr/bin/env bash
#SBATCH -J hg002_prepare_eval
#SBATCH -o 04_prepare_eval.log
#SBATCH -p rl8
#SBATCH -c 64
#SBATCH -n 1
#SBATCH -N 1
set -euo pipefail
source config.sh
module load ${_MINIMAP2}
set -x
THREADS=${SLURM_CPUS_PER_TASK:-64}

# Resolvability
python3 "$SCRIPTS/resolvable.py" --threads "$THREADS" --near 2 --med 10 \
  "$REF" "$READS_CONFIDENT_LOCUS" "$D_EVAL/resolvability_strata.tsv"

# Hapmer strata
awk 'BEGIN{OFS="\t"}
 NR==FNR {map_hap[$1]=$2; next}
 !($1 in map_hap) {next}
 {h=$5+$6
  if(h==0) s="no_hapmer"
  else if($2!="p" && $2!="m") s="ambiguous"
  else if(($2=="p"?"pat":"mat")!=map_hap[$1]) s="disagree"
  else s=(h<=40?"vhard":(h<=100?"hard":(h<=350?"med":"easy")))
  print $1,s}' "$READS_CONFIDENT_LOCUS" "$D_TRIO/$READS_BASE.triobin.tsv" > "$D_EVAL/hapmer_strata.tsv"

# Variant strata
python3 "$SCRIPTS/variant_strata.py" \
  "$D_REF/variants.tsv" "$READS_CONFIDENT_LOCUS" "$D_EVAL/variant_strata.tsv"
