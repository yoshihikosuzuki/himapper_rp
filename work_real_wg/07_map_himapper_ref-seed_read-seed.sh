#!/usr/bin/env bash
#SBATCH -J hg002_wg_himapper
#SBATCH -o 07_map_himapper_ref-seed_read-seed.%j.log
#SBATCH -p rl8
#SBATCH -c 16
#SBATCH -n 1
#SBATCH -N 1
set -euo pipefail
source config.sh
module load ${_DALIGNER} ${_PYTHON}
set -x
THREADS=${SLURM_CPUS_PER_TASK:-16}
HIM="$REPO/himapper"
RUN_DIR="$D_HIMAPPER/ref-seed_read-seed"
HIMAPPER_TMP="$WG_TMP_ROOT/himapper"
mkdir -p "$RUN_DIR" "$HIMAPPER_TMP"
cd "$RUN_DIR"

sha256sum "$HIM" "$DAZZ_DB_ROOT/DBshow" \
  "$SCRIPTS/lashow_to_placement.py"
ln -sf "$D_HIMAPPER"/{ref.dam,.ref.bps,.ref.hdr,.ref.idx} .
ln -sf "$D_HIMAPPER"/{reads.dam,.reads.bps,.reads.hdr,.reads.idx} .
NB=$(awk -F'=' '/blocks/{gsub(/ /,"",$2);print $2}' reads.dam)
export CN_UNIQ_W CN_TIER_SRC CN_HASH_TIE CN_DUMP

MAP_START=${HIMAPPER_MAP_START:-map}
printf 'HIMAPPER_MAP_START=%s\n' "$MAP_START"

if [[ "$MAP_START" == map ]]; then
  if [[ ! "$THREADS" =~ ^[0-9]{1,3}$ ]]; then
    echo "Integrated himapper requires an integer thread count, got: $THREADS" >&2
    exit 1
  fi
  THREADS=$(( 10#$THREADS ))
  if (( THREADS < 1 || THREADS > 128 )); then
    echo "Integrated himapper requires 1..128 total threads, got: $THREADS" >&2
    exit 1
  fi

  SEED_MODE=${HIMAPPER_SEED_MODE:-pipeline}
  SEED_THREAD_OVERRIDE=${HIMAPPER_SEED_THREADS:-}
  case "$SEED_MODE" in
    serial)
      if [[ -n "$SEED_THREAD_OVERRIDE" ]]; then
        if [[ ! "$SEED_THREAD_OVERRIDE" =~ ^[0-9]{1,3}$ ]] ||
           (( 10#$SEED_THREAD_OVERRIDE != THREADS )); then
          echo "Serial mode uses the full T/T budget; HIMAPPER_SEED_THREADS must be unset or $THREADS" >&2
          exit 1
        fi
      fi
      SEED_THREADS=$THREADS
      MAP_THREADS=$THREADS
      ;;
    pipeline)
      if (( THREADS < 2 )); then
        echo "Pipeline mode requires at least 2 total threads" >&2
        exit 1
      fi
      if [[ -n "$SEED_THREAD_OVERRIDE" ]]; then
        if [[ ! "$SEED_THREAD_OVERRIDE" =~ ^[0-9]{1,3}$ ]]; then
          echo "HIMAPPER_SEED_THREADS must be a positive integer" >&2
          exit 1
        fi
        SEED_THREADS=$(( 10#$SEED_THREAD_OVERRIDE ))
      else
        # Nearest integer to 3/8 of the total: T=16 gives the measured 6/10 split.
        SEED_THREADS=$(( (3 * THREADS + 4) / 8 ))
        (( SEED_THREADS < 1 )) && SEED_THREADS=1
        (( SEED_THREADS >= THREADS )) && SEED_THREADS=$(( THREADS - 1 ))
      fi
      if (( SEED_THREADS < 1 || SEED_THREADS >= THREADS )); then
        echo "Pipeline split must satisfy 1 <= seed threads < total threads; got $SEED_THREADS/$THREADS" >&2
        exit 1
      fi
      MAP_THREADS=$(( THREADS - SEED_THREADS ))
      ;;
    *)
      echo "Invalid HIMAPPER_SEED_MODE=$SEED_MODE; use serial or pipeline" >&2
      exit 1
      ;;
  esac

  SEED_ARGS=(
    --read-seed-profiles "$D_HIMAPPER/reads.ref.prof" "$D_HIMAPPER/reads.prof"
    --seed-mode "$SEED_MODE"
    --seed-threads "$SEED_THREADS"
    --seed-window "$HIMAPPER_WS"
    --seed-spacing "$HIMAPPER_WS"
  )
  if [[ -n "${HIMAPPER_SAVE_READ_KPOS:-}" ]]; then
    SEED_ARGS+=(--save-read-kpos "$HIMAPPER_SAVE_READ_KPOS")
  fi
  printf 'HIMAPPER_SEED_MODE=%s total_threads=%d seed_threads=%d map_threads=%d save_read_kpos=%s\n' \
    "$SEED_MODE" "$THREADS" "$SEED_THREADS" "$MAP_THREADS" \
    "${HIMAPPER_SAVE_READ_KPOS:-off}"
fi

case "$MAP_START" in
  map)
    run_timed "$RUN_DIR/seed-map.time" \
      "$HIM" -k"$KMER" -1 -v -T"$THREADS" -P"$HIMAPPER_TMP" \
      "${SEED_ARGS[@]}" \
      ref.dam "$D_HIMAPPER/.ref.minimax.kpos" \
      reads.dam generate
    LAcat "reads.@1-${NB}.ref.las" > "$RUN_DIR/reads.ref.las"
    ;;
  placement)
    if [[ ! -s "$RUN_DIR/reads.ref.las" ]]; then
      echo "Cannot resume placement: missing $RUN_DIR/reads.ref.las" >&2
      exit 1
    fi
    ;;
  *)
    echo "Invalid HIMAPPER_MAP_START=$MAP_START; use map or placement" >&2
    exit 1
    ;;
esac

REF_LAYOUT="$RUN_DIR/ref.dam.layout"
MAPPING_TMP="$RUN_DIR/mapping.tsv.${SLURM_JOB_ID:?}.tmp"
"$DAZZ_DB_ROOT/DBshow" -n ref.dam > "$REF_LAYOUT"
trap 'rm -f "$MAPPING_TMP"' EXIT
LAshow reads ref.dam "$RUN_DIR/reads.ref.las" 2>/dev/null \
  | python3 "$SCRIPTS/lashow_to_placement.py" \
      --ref-layout "$REF_LAYOUT" \
      "$D_HIMAPPER/reads.names" "$D_HIMAPPER/ref.names" \
  > "$MAPPING_TMP"
mv -f "$MAPPING_TMP" "$RUN_DIR/mapping.tsv"
trap - EXIT
