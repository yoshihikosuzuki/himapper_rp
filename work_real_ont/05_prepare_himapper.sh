#!/usr/bin/env bash
#SBATCH -J hg002_ont_himapper_prep
#SBATCH -o 05_prepare_himapper.log
#SBATCH -p rl8
#SBATCH -c 16
#SBATCH -n 1
#SBATCH -N 1
set -euo pipefail
source config.sh
for tool in "$DAZZ_DB_ROOT/fasta2DAM" "$DAZZ_DB_ROOT/DBshow" \
            "$DAZZ_DB_ROOT/DBsplit" \
            "$FASTK_ROOT/FastK" "$FASTK_ROOT/Fastrm" \
            "$REPO/fastkdual/FastKdual"; do
  if [[ ! -x "$tool" ]]; then
    echo "Missing local tool: $tool; run 'make local-tools' in $REPO" >&2
    exit 1
  fi
done
set -x
THREADS=${SLURM_CPUS_PER_TASK:-16}
FASTK_TMP="$ONT_TMP_ROOT/fastk"
mkdir -p "$FASTK_TMP"

cd "$D_HIMAPPER"
sha256sum "$DAZZ_DB_ROOT/fasta2DAM" "$DAZZ_DB_ROOT/DBshow" \
  "$DAZZ_DB_ROOT/DBsplit" \
  "$FASTK_ROOT/FastK" "$FASTK_ROOT/Fastrm" \
  "$REPO/fastkdual/FastKdual" "$REPO/minimaxmer"

"$DAZZ_DB_ROOT/fasta2DAM" ref "$REF"
"$DAZZ_DB_ROOT/fasta2DAM" reads "$READS"
"$DAZZ_DB_ROOT/DBsplit" -s"$DB_BLOCK_SIZE" -f ref.dam
"$DAZZ_DB_ROOT/DBsplit" -s"$DB_BLOCK_SIZE" -f reads.dam

ln -sf "$READ_NAMES" reads.names
cut -f1 "${REF}.fai" > ref.names

run_timed "$D_HIMAPPER/ref-profile.time" \
  "$FASTK_ROOT/FastK" -k"$KMER" -v -T"$THREADS" -t1 -p -P"$FASTK_TMP" ref.dam
run_timed "$D_HIMAPPER/reads-dual-profile.time" \
  "$REPO/fastkdual/FastKdual" -k"$KMER" -v -T"$THREADS" -t1 \
  -p:ref -Nreads -Rreads.ref -P"$FASTK_TMP" reads.dam
run_timed "$D_HIMAPPER/ref-vs-reads-profile.time" \
  "$FASTK_ROOT/FastK" -k"$KMER" -v -T"$THREADS" -t1 -p:reads -Nref.reads \
  -P"$FASTK_TMP" ref.dam

run_timed "$D_HIMAPPER/ref-seeds.time" \
  "$REPO/minimaxmer" -T"$THREADS" -k"$KMER" \
  -w"$HIMAPPER_WS" -s"$HIMAPPER_WS" \
  ref.dam ref.prof ref.reads.prof "$D_HIMAPPER/.ref.minimax.kpos"
run_timed "$D_HIMAPPER/read-seeds.time" \
  "$REPO/minimaxmer" -T"$THREADS" -k"$KMER" \
  -w"$HIMAPPER_WS" -s"$HIMAPPER_WS" \
  reads.dam reads.ref.prof reads.prof "$D_HIMAPPER/.reads.minimax.kpos"
