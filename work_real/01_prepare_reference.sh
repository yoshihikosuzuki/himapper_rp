#!/usr/bin/env bash
#SBATCH -J hg002_ref
#SBATCH -o 01_prepare_reference.log
#SBATCH -p rl8
#SBATCH -c 16
#SBATCH -n 1
#SBATCH -N 1
set -euo pipefail
source config.sh
# Keep the exact v1.1 FASTA used by historical mappings and indexes.
if [[ -L "$REF_ALL" || -L "$REF" ]]; then
  printf 'HG002 reference FASTAs are centralized; use the existing links instead of rebuilding in place.\n' >&2
  exit 1
fi
module load ${_SAMTOOLS} ${_MINIMAP2}
set -x
THREADS=${SLURM_CPUS_PER_TASK:-16}

# Build diploid reference
( zcat "$D_REF/hg002.mat.fa.gz"
  zcat "$D_REF/hg002.pat.fa.gz" ) > "$REF_ALL"
samtools faidx "$REF_ALL"

# Extract chr6 haplotypes
samtools faidx "$REF_ALL" "$CHR6_MAT" "$CHR6_PAT" > "$REF"
samtools faidx "$REF"

# Align chr6 haplotypes
minimap2 -c --cs=long -x asm5 -t"$THREADS" \
  <(samtools faidx "$REF" "$CHR6_PAT") \
  <(samtools faidx "$REF" "$CHR6_MAT") \
  > "$D_REF/mat_vs_pat.paf"
minimap2 -c --cs=long -x asm5 -t"$THREADS" \
  <(samtools faidx "$REF" "$CHR6_MAT") \
  <(samtools faidx "$REF" "$CHR6_PAT") \
  > "$D_REF/pat_vs_mat.paf"

# Extract MAT/PAT variants
python3 "$SCRIPTS/gen_variants.py" "$D_REF/variants.tsv" "$D_REF/mat_vs_pat.paf" "$D_REF/pat_vs_mat.paf"
