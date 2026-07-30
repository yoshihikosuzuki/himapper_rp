#!/usr/bin/env bash
#SBATCH -J ref_seed_secondary
#SBATCH -o 09_analyze_secondary.log
#SBATCH -p rl8
#SBATCH -c 64
#SBATCH -n 1
#SBATCH -N 1
# Separate candidate/secondary analysis for ref-seed/read-seed at W=11.5.
# This is not part of the primary-only benchmark.
set -euo pipefail
source config.sh
module load ${_DALIGNER}
set -x

THREADS=${SLURM_CPUS_PER_TASK:-64}
HIM="$REPO/himapper"
RUN_DIR="$D_HIMAPPER/ref-seed_read-seed-secondary"
WORK_DIR="$RUN_DIR/map-work"
mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

ln -sf "$D_HIMAPPER"/{ref.dam,.ref.bps,.ref.hdr,.ref.idx} .
ln -sf "$D_HIMAPPER"/{reads.dam,.reads.bps,.reads.hdr,.reads.idx} .

NB=$(awk -F'=' '/blocks/{gsub(/ /,"",$2);print $2}' reads.dam)
REF_LAYOUT="$RUN_DIR/ref.dam.layout"
"$DAZZ_DB_ROOT/DBshow" -n ref.dam > "$REF_LAYOUT"
export CN_UNIQ_W=11.5 CN_TIER_SRC=1 CN_HASH_TIE=1 CN_DUMP=0
"$HIM" -n0.7 -T"$THREADS" \
  ref.dam "$D_HIMAPPER/.ref.minimax.kpos" \
  reads.dam "$D_HIMAPPER/.reads.minimax.kpos"
LAcat "reads.@1-${NB}.ref.las" > "$RUN_DIR/reads.ref.las"
LAshow reads ref.dam "$RUN_DIR/reads.ref.las" 2>/dev/null \
  | python3 "$SCRIPTS/lashow_to_placement.py" --all \
      --ref-layout "$REF_LAYOUT" \
      "$D_HIMAPPER/reads.names" "$D_HIMAPPER/ref.names" \
  > "$RUN_DIR/mapping.all-candidates.tsv"
# Headerless candidate-chain-v1 TSV: read, hap, scaffold, 0-based half-open
# reference start/end, strand, and reference bounding span (internal gaps
# included).  Each row is one bounded START/NEXT chain in LAshow input order.
# Unmapped reads are omitted and the default-selected chain is not marked.
# This is for candidate diagnostics only, not primary eval/coverage.
awk -F'\t' '
  { c[$1]++ }
  END {
    for (r in c) {
      n = c[r]
      if (n >= 6) d[6]++
      else d[n]++
    }
    for (n = 1; n <= 5; n++)
      printf "  %d cand: %d reads\n", n, d[n] + 0
    printf "  6+ cand: %d reads\n", d[6] + 0
  }
' \
  "$RUN_DIR/mapping.all-candidates.tsv" >&2
