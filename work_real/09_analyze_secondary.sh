#!/usr/bin/env bash
#SBATCH -J ref_seed_secondary
#SBATCH -o 09_analyze_secondary.log
#SBATCH -p rl8
#SBATCH -c 64
#SBATCH -n 1
#SBATCH -N 1
# Separate candidate/secondary analysis for ref-seed/read-all at W=11.5.
# This is not part of the primary-only benchmark.
set -euo pipefail
source config.sh
module load ${_DALIGNER}
set -x

THREADS=${SLURM_CPUS_PER_TASK:-64}
HIM="$REPO/himapper"
HPC="$REPO/HPC.himapper"
RUN_DIR="$D_HIMAPPER/ref-seed_read-all-secondary"
WORK_DIR="$RUN_DIR/map-work"
mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

ln -sf "$D_HIMAPPER"/{ref.dam,.ref.bps,.ref.hdr,.ref.idx} .
ln -sf "$D_HIMAPPER"/{reads.dam,.reads.bps,.reads.hdr,.reads.idx} .

NB=$(awk -F'=' '/blocks/{gsub(/ /,"",$2);print $2}' reads.dam)
export CN_UNIQ_W=11.5 CN_TIER_SRC=0 CN_HASH_TIE=1 CN_DUMP=0
COMMANDS="$RUN_DIR/commands.sh"
{
  echo '#!/usr/bin/env bash'
  echo 'set -euo pipefail'
  "$HPC" -n0.7 -T"$THREADS" -K"$D_HIMAPPER/.ref.#.minimax.kpos" ref.dam reads 2>/dev/null \
    | grep '^himapper' | sed "s|^himapper|$HIM|"
} > "$COMMANDS"
[[ $(tail -n +3 "$COMMANDS" | wc -l) -eq "$NB" ]] || {
  echo "generated command count does not match blocks=$NB" >&2
  exit 1
}
rm -f reads.*.ref.las
bash "$COMMANDS"
: > "$RUN_DIR/mapping.all-candidates.tsv"
for k in $(seq 1 "$NB"); do
  LAshow reads ref.dam "reads.$k.ref.las" 2>/dev/null \
    | python3 "$SCRIPTS/lashow_to_placement.py" --all "$D_HIMAPPER/reads.names" \
        "$D_HIMAPPER/ref.names" >> "$RUN_DIR/mapping.all-candidates.tsv"
done
# Keep the LAS files: unlike mapping.all-candidates.tsv, they retain the
# chain flags needed to distinguish candidates from split alignment pieces.
awk -F'\t' '{c[$1]++} END{for(r in c)d[c[r]]++; for(n=1;n<=5;n++)printf "  %d cand: %d reads\n",n,d[n]+0; s=0;for(n=6;n<=999;n++)s+=d[n];printf "  6+ cand: %d reads\n",s}' \
  "$RUN_DIR/mapping.all-candidates.tsv" >&2
