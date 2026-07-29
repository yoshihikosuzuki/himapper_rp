#!/usr/bin/env bash
#SBATCH -J hg002_ont_reads
#SBATCH -o 02_prepare_reads.log
#SBATCH -p rl8
#SBATCH -c 192
#SBATCH -n 1
#SBATCH -N 1
set -euo pipefail
source config.sh
module load ${_SEQKIT} ${_WINNOWMAP}
set -x
THREADS=${SLURM_CPUS_PER_TASK:-192}
SORT_TMP="$ONT_TMP_ROOT/truth-sort"
mkdir -p "$SORT_TMP"

# Assign truth against the complete diploid assembly before selecting chr6.
run_timed "$D_READS/truth-map.time" \
  winnowmap -I8G -W "$DIP_REPETITIVE" -x map-ont \
    -t"$THREADS" --secondary=no -c "$REF_ALL" "${ONT_FASTQS[@]}" \
  > "$TRUTH_SOURCE_PAF"

# Resolve split records by the largest aligned query span, then retain reads
# whose best placement is on maternal or paternal chr6.
awk -v min_len="$MIN_READ_LENGTH" -v min_cov="$MIN_QUERY_COVERAGE" \
  'BEGIN{OFS="\t"}
   {cov=($4-$3)/$2
    if($2>min_len && cov>=min_cov){
      hap=($6 ~ /_MATERNAL$/ ? "mat" :
           ($6 ~ /_PATERNAL$/ ? "pat" : "NA"))
      print $1,hap,$6,$8,$9,($4-$3),$12
    }}' "$TRUTH_SOURCE_PAF" \
  | LC_ALL=C sort -T "$SORT_TMP" -k1,1 -k6,6nr \
  | awk '!seen[$1]++' \
  | awk -v m="$CHR6_MAT" -v p="$CHR6_PAT" '$3==m||$3==p' \
  > "$READS_MAPPING_TRUTH"

awk -F'\t' -v min_mapq="$MIN_TRUTH_MAPQ" '$7>=min_mapq' \
  "$READS_MAPPING_TRUTH" > "$READS_CONFIDENT_MAPPING_TRUTH"

seqkit grep --quiet -j "$THREADS" \
  -f <(cut -f1 "$READS_MAPPING_TRUTH") "${ONT_FASTQS[@]}" \
  | seqkit fq2fa --quiet -j "$THREADS" -w0 > "$READS"

seqkit seq --quiet -j "$THREADS" -n -i "$READS" > "$READ_NAMES"
