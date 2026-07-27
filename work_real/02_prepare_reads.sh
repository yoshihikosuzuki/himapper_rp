#!/usr/bin/env bash
#SBATCH -J hg002_reads
#SBATCH -o 02_prepare_reads.log
#SBATCH -p rl8
#SBATCH -c 192
#SBATCH -n 1
#SBATCH -N 1
set -euo pipefail
source config.sh
module load ${_MERYL} ${_WINNOWMAP} ${_SEQKIT}
set -x
THREADS=${SLURM_CPUS_PER_TASK:-192}

DIP_MERYL="$D_READS/dip.meryl"
DIP_REP="$D_READS/dip.repetitive.txt"
READS_PAF="$D_READS/reads.map.paf"

meryl count k=15 output "$DIP_MERYL" "$REF_ALL"
meryl print greater-than distinct=0.9998 "$DIP_MERYL" > "$DIP_REP"

winnowmap -I8G -W "$DIP_REP" -x map-pb -t"$THREADS" --secondary=no -c \
  "$REF_ALL" "$READS_ALL" > "$READS_PAF"

# Select one chr6 mapping-truth locus per read
awk 'BEGIN{OFS="\t"}
     { cov=($4-$3)/$2;
       if($2>=10000 && cov>=0.90){
         hap=($6 ~ /_MATERNAL$/ ? "mat" : "pat");
         print $1, hap, $6, $8, $9, ($4-$3), $12 } }' \
  "$READS_PAF" \
  | LC_ALL=C sort -k1,1 -k6,6nr \
  | awk '!seen[$1]++' \
  | awk -v m="$CHR6_MAT" -v p="$CHR6_PAT" '$3==m||$3==p' \
  > "$READS_MAPPING_TRUTH"

# High-confidence mapping truth for stratification and accuracy evaluation
awk -F'\t' '$7>=60' "$READS_MAPPING_TRUTH" \
  > "$READS_CONFIDENT_MAPPING_TRUTH"

# Extract selected reads
seqkit grep --quiet -j "$THREADS" \
  -f <(cut -f1 "$READS_MAPPING_TRUTH") \
  "$READS_ALL" \
  | seqkit fq2fa --quiet -j "$THREADS" -w0 > "$READS"

seqkit seq --quiet -j "$THREADS" -n -i "$READS" > "$READ_NAMES"
