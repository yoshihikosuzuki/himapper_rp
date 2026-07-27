#!/usr/bin/env bash
#SBATCH -J hg002_prepare_himapper
#SBATCH -o 05_prepare_himapper.log
#SBATCH -p rl8
#SBATCH -c 16
#SBATCH -n 1
#SBATCH -N 1
set -euo pipefail
source config.sh
module load ${_DAZZ_DB} ${_FASTK}
set -x
THREADS=${SLURM_CPUS_PER_TASK:-16}
# Reuse existing DAM/profiles while re-splitting the DBs and regenerating KPOS.
SEEDS_ONLY=${SEEDS_ONLY:-0}

cd "$D_HIMAPPER"

if [[ "$SEEDS_ONLY" != 1 ]]; then
  # Build DAMs in a fresh preparation directory.
  fasta2DAM ref "$REF"
  fasta2DAM reads "$READS"
fi

# KPOS runs and the 10,000-hit threshold follow these DB block boundaries.
DBsplit -s200 -f ref.dam
DBsplit -s200 -f reads.dam

if [[ "$SEEDS_ONLY" != 1 ]]; then
  ln -sf "$READ_NAMES" reads.names
  cut -f1 "${REF}.fai" > ref.names

  run_timed "$D_HIMAPPER/ref-profile.time" \
    FastK -k"$KMER" -v -T"$THREADS" -t1 -p -P"$TMPDIR" ref.dam
  run_timed "$D_HIMAPPER/reads-profile.time" \
    FastK -k"$KMER" -v -T"$THREADS" -t1 -p -P"$TMPDIR" reads.dam
  run_timed "$D_HIMAPPER/reads-vs-ref-profile.time" \
    FastK -k"$KMER" -v -T"$THREADS" -t1 -p:ref -Nreads.ref -P"$TMPDIR" reads.dam
  run_timed "$D_HIMAPPER/ref-vs-reads-profile.time" \
    FastK -k"$KMER" -v -T"$THREADS" -t1 -p:reads -Nref.reads -P"$TMPDIR" ref.dam
fi

# Generate seeds
run_timed "$D_HIMAPPER/ref-seeds.time" \
  "$REPO/minimaxmer" -T"$THREADS" ref.dam ref.prof ref.reads.prof \
  "$D_HIMAPPER/.ref.minimax.kpos"
run_timed "$D_HIMAPPER/read-seeds.time" \
  "$REPO/minimaxmer" -T"$THREADS" reads.dam reads.ref.prof reads.prof \
  "$D_HIMAPPER/.reads.minimax.kpos"
