#!/usr/bin/env bash
#SBATCH -J hg002_wg_himapper_prep
#SBATCH -o 05_prepare_himapper.%j.log
#SBATCH -p rl8
#SBATCH -c 16
#SBATCH -n 1
#SBATCH -N 1
set -euo pipefail
source config.sh
for tool in "$DAZZ_DB_ROOT/fasta2DAM" "$DAZZ_DB_ROOT/DBshow" \
            "$DAZZ_DB_ROOT/DBsplit" \
            "$DAZZ_DB_ROOT/DBstats" \
            "$FASTK_ROOT/FastK" "$FASTK_ROOT/Fastrm"; do
  if [[ ! -x "$tool" ]]; then
    echo "Missing local tool: $tool; run 'make local-tools' in $REPO" >&2
    exit 1
  fi
done
set -x
THREADS=${SLURM_CPUS_PER_TASK:-16}
FASTK_TMP="$WG_TMP_ROOT/fastk"
mkdir -p "$FASTK_TMP"

PREP_START=${HIMAPPER_PREP_START:-dam}
case "$PREP_START" in
  dam)      PREP_START_RANK=0 ;;
  split)    PREP_START_RANK=1 ;;
  profiles) PREP_START_RANK=2 ;;
  seeds)    PREP_START_RANK=3 ;;
  *)
    echo "Invalid HIMAPPER_PREP_START=$PREP_START; use dam, split, profiles, or seeds" >&2
    exit 1
    ;;
esac

cd "$D_HIMAPPER"
sha256sum "$DAZZ_DB_ROOT/fasta2DAM" "$DAZZ_DB_ROOT/DBshow" \
  "$DAZZ_DB_ROOT/DBsplit" \
  "$DAZZ_DB_ROOT/DBstats" \
  "$FASTK_ROOT/FastK" "$FASTK_ROOT/Fastrm" "$REPO/minimaxmer"
printf 'HIMAPPER_PREP_START=%s HIMAPPER_WS=%s DB_BLOCK_SIZE=%s\n' \
  "$PREP_START" "$HIMAPPER_WS" "$DB_BLOCK_SIZE"

verify_dam()
{
  local root=$1
  local input=$2
  local expected_source
  local actual_source
  local nfiles
  local path

  for path in "${root}.dam" ".${root}.idx" ".${root}.bps" ".${root}.hdr"; do
    if [[ ! -s "$path" ]]; then
      echo "Cannot resume: missing or empty $D_HIMAPPER/$path" >&2
      exit 1
    fi
  done
  if [[ ! -s "$input" ]]; then
    echo "Cannot resume: missing or empty input $input" >&2
    exit 1
  fi
  if [[ "$input" -nt "${root}.dam" ]]; then
    echo "Cannot resume: input $input is newer than ${root}.dam" >&2
    exit 1
  fi

  expected_source=$(basename "$input")
  expected_source=${expected_source%.gz}
  expected_source=${expected_source%.fasta}
  expected_source=${expected_source%.fastq}
  expected_source=${expected_source%.fa}
  expected_source=${expected_source%.fq}
  nfiles=$(awk '$1 == "files" {print $3; exit}' "${root}.dam")
  actual_source=$(awk 'NR == 2 {print $2; exit}' "${root}.dam")
  if [[ "$nfiles" != 1 || "$actual_source" != "$expected_source" ]]; then
    echo "Cannot resume: ${root}.dam does not match input $input" >&2
    exit 1
  fi
}

create_dam()
{
  local root=$1
  local input=$2
  local path

  for path in "${root}.dam" ".${root}.idx" ".${root}.bps" ".${root}.hdr"; do
    if [[ -e "$path" ]]; then
      echo "Refusing to append to existing $D_HIMAPPER/$path;" \
           "set HIMAPPER_PREP_START=split to reuse a completed DAM" >&2
      exit 1
    fi
  done
  "$DAZZ_DB_ROOT/fasta2DAM" "$root" "$input"
  verify_dam "$root" "$input"
}

verify_split()
{
  local root=$1
  local expected_size=$((DB_BLOCK_SIZE * 1000000))
  local blocks
  local size

  blocks=$(awk '$1 == "blocks" {print $3; exit}' "${root}.dam")
  size=$(awk '$1 == "size" {print $3; exit}' "${root}.dam")
  if [[ ! "$blocks" =~ ^[1-9][0-9]*$ || "$size" != "$expected_size" ]]; then
    echo "Invalid split in ${root}.dam: blocks=$blocks size=$size;" \
         "expected a positive block count and size=$expected_size" >&2
    exit 1
  fi
}

verify_profiles()
{
  local root
  local part
  local parts

  shopt -s nullglob
  for root in ref reads reads.ref ref.reads; do
    if [[ ! -s "${root}.prof" ]]; then
      echo "Missing FastK profile: $D_HIMAPPER/${root}.prof" >&2
      exit 1
    fi
    parts=(".${root}.pidx."*)
    if (( ${#parts[@]} == 0 )); then
      echo "Missing FastK profile index parts for $root" >&2
      exit 1
    fi
    for part in "${parts[@]}"; do
      [[ -s "$part" ]] || {
        echo "Empty FastK profile index part: $part" >&2
        exit 1
      }
    done
    parts=(".${root}.prof."*)
    if (( ${#parts[@]} == 0 )); then
      echo "Missing FastK profile data parts for $root" >&2
      exit 1
    fi
    for part in "${parts[@]}"; do
      [[ -s "$part" ]] || {
        echo "Empty FastK profile data part: $part" >&2
        exit 1
      }
    done
  done
  shopt -u nullglob
}

if (( PREP_START_RANK == 0 )); then
  create_dam ref "$REF"
  create_dam reads "$READS"
else
  verify_dam ref "$REF"
  verify_dam reads "$READS"
fi

if (( PREP_START_RANK <= 1 )); then
  "$DAZZ_DB_ROOT/DBstats" ref.dam \
    > "ref.pre-split.${SLURM_JOB_ID}.dbstats.txt"
  "$DAZZ_DB_ROOT/DBstats" reads.dam \
    > "reads.pre-split.${SLURM_JOB_ID}.dbstats.txt"
  "$DAZZ_DB_ROOT/DBsplit" -s"$DB_BLOCK_SIZE" -f ref.dam
  "$DAZZ_DB_ROOT/DBsplit" -s"$DB_BLOCK_SIZE" -f reads.dam
fi

verify_split ref
verify_split reads
"$DAZZ_DB_ROOT/DBstats" ref.dam \
  > "ref.post-split.${SLURM_JOB_ID}.dbstats.txt"
"$DAZZ_DB_ROOT/DBstats" reads.dam \
  > "reads.post-split.${SLURM_JOB_ID}.dbstats.txt"

ln -sf "$READ_NAMES" reads.names
cut -f1 "${REF}.fai" > ref.names

if (( PREP_START_RANK <= 2 )); then
  run_timed "$D_HIMAPPER/ref-profile.time" \
    "$FASTK_ROOT/FastK" -k"$KMER" -v -T"$THREADS" -t1 -p -P"$FASTK_TMP" ref.dam
  run_timed "$D_HIMAPPER/reads-profile.time" \
    "$FASTK_ROOT/FastK" -k"$KMER" -v -T"$THREADS" -t1 -p -P"$FASTK_TMP" reads.dam
  run_timed "$D_HIMAPPER/reads-vs-ref-profile.time" \
    "$FASTK_ROOT/FastK" -k"$KMER" -v -T"$THREADS" -t1 -p:ref -Nreads.ref \
    -P"$FASTK_TMP" reads.dam
  run_timed "$D_HIMAPPER/ref-vs-reads-profile.time" \
    "$FASTK_ROOT/FastK" -k"$KMER" -v -T"$THREADS" -t1 -p:reads -Nref.reads \
    -P"$FASTK_TMP" ref.dam
fi

verify_profiles

run_timed "$D_HIMAPPER/ref-seeds.time" \
  "$REPO/minimaxmer" -T"$THREADS" -k"$KMER" \
  -w"$HIMAPPER_WS" -s"$HIMAPPER_WS" \
  ref.dam ref.prof ref.reads.prof "$D_HIMAPPER/.ref.minimax.kpos"
run_timed "$D_HIMAPPER/read-seeds.time" \
  "$REPO/minimaxmer" -T"$THREADS" -k"$KMER" \
  -w"$HIMAPPER_WS" -s"$HIMAPPER_WS" \
  reads.dam reads.ref.prof reads.prof "$D_HIMAPPER/.reads.minimax.kpos"
