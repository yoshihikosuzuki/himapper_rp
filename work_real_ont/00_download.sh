#!/usr/bin/env bash
#SBATCH -J hg002_ont_download
#SBATCH -o 00_download.log
#SBATCH -p rl8
#SBATCH -c 4
#SBATCH -n 1
#SBATCH -N 1
set -euo pipefail
source config.sh
module load ${_SEQKIT}
set -x

for i in "${!ONT_FASTQ_NAMES[@]}"; do
  wget -c -P "$D_RAW" "$ONT_FASTQ_BASE/${ONT_FASTQ_NAMES[$i]}"
  test "$(stat -c%s "${ONT_FASTQS[$i]}")" -eq "${ONT_FASTQ_SIZES[$i]}"
done

{
  printf 'file\tsampled_reads\tmin_length\tmax_length\tle_100kb\tgt_100kb\n'
  for fastq in "${ONT_FASTQS[@]}"; do
    seqkit head --quiet -n 100 "$fastq" \
      | seqkit fx2tab --quiet -n -l \
      | awk -v file="$(basename "$fastq")" '
          NR==1 {min=$2; max=$2}
          {
            if($2<min) min=$2
            if($2>max) max=$2
            if($2<=100000) le++
            else gt++
          }
          END {
            print file,NR,min,max,le+0,gt+0
          }' OFS='\t'
  done
} > "$D_RAW/first100_read_lengths.tsv"
