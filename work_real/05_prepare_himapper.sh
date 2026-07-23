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

cd "$D_HIMAPPER"

# Build DAMs and profiles
fasta2DAM ref "$REF"
run_timed "$D_HIMAPPER/ref-profile.time" \
  FastK -k"$KMER" -v -T"$THREADS" -t1 -p -P"$TMPDIR" ref.dam

fasta2DAM reads "$READS"
DBsplit -s200 -f reads.dam

cp "$READ_NAMES" reads.names
cut -f1 "${REF}.fai" > ref.names

run_timed "$D_HIMAPPER/reads-profile.time" \
  FastK -k"$KMER" -v -T"$THREADS" -t1 -p -P"$TMPDIR" reads.dam
run_timed "$D_HIMAPPER/reads-vs-ref-profile.time" \
  FastK -k"$KMER" -v -T"$THREADS" -t1 -p:ref -Nreads.ref -P"$TMPDIR" reads.dam
run_timed "$D_HIMAPPER/ref-vs-reads-profile.time" \
  FastK -k"$KMER" -v -T"$THREADS" -t1 -p:reads -Nref.reads -P"$TMPDIR" ref.dam

# Generate seeds
run_timed "$D_HIMAPPER/ref-seeds.time" \
  "$REPO/minimaxmer" ref.dam ref.prof ref.reads.prof "$D_HIMAPPER/ref.minimax.kpos"
run_timed "$D_HIMAPPER/read-seeds.time" \
  "$REPO/minimaxmer" reads.dam reads.ref.prof reads.prof "$D_HIMAPPER/reads.minimax.kpos"
