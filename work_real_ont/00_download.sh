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
