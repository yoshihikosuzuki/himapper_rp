#!/usr/bin/env bash
#SBATCH -J hg002_wg_reads
#SBATCH -o 02_prepare_reads.log
#SBATCH -p rl8
#SBATCH -c 64
#SBATCH -n 1
#SBATCH -N 1
set -euo pipefail
source config.sh
module load ${_SEQKIT}
set -x
THREADS=${SLURM_CPUS_PER_TASK:-64}
SORT_TMP="$WG_TMP_ROOT/truth-sort"
mkdir -p "$SORT_TMP"

# Select one whole-genome mapping-truth locus per read. As in work_real,
# split records are resolved by the largest aligned query span.
awk -v min_len="$MIN_READ_LENGTH" -v min_cov="$MIN_QUERY_COVERAGE" \
  'BEGIN{OFS="\t"}
   {cov=($4-$3)/$2
    if($2>=min_len && cov>=min_cov){
      hap=($6 ~ /_MATERNAL$/ ? "mat" : "pat")
      print $1,hap,$6,$8,$9,($4-$3),$12
    }}' "$TRUTH_SOURCE_PAF" \
  | LC_ALL=C sort -T "$SORT_TMP" -k1,1 -k6,6nr \
  | awk '!seen[$1]++' \
  > "$READS_MAPPING_TRUTH"

awk -F'\t' -v min_mapq="$MIN_TRUTH_MAPQ" '$7>=min_mapq' \
  "$READS_MAPPING_TRUTH" > "$READS_CONFIDENT_MAPPING_TRUTH"

seqkit grep --quiet -j "$THREADS" \
  -f <(cut -f1 "$READS_MAPPING_TRUTH") "$READS_ALL" \
  | seqkit fq2fa --quiet -j "$THREADS" -w0 > "$READS"

seqkit seq --quiet -j "$THREADS" -n -i "$READS" > "$READ_NAMES"
