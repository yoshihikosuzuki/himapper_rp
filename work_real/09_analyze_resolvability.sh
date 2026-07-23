#!/usr/bin/env bash
#SBATCH -J reslvl
#SBATCH -o 09_analyze_resolvability.log
#SBATCH -p rl8
#SBATCH -c 1
#SBATCH -n 1
#SBATCH -N 1
# Graded resolvability (L0..L3) x agreement for himapper (W=11.5) / mm2 / wm2.
set -euo pipefail
source config.sh
set -x
python3 "$SCRIPTS/eval.py" --loci "$READS_CONFIDENT_LOCUS" \
    --resolvable-mode kmer --resolvable "$D_EVAL/resolvability_strata.tsv" \
    --strata "$D_EVAL/resolvability_strata.tsv" --ov 0.5 \
    --liftover-paf "$D_REF/mat_vs_pat.paf" \
    --liftover-paf "$D_REF/pat_vs_mat.paf" \
    ref-seed_read-all="$D_HIMAPPER/ref-seed_read-all/placement.tsv" \
    mm2="$D_MM2/placement.tsv" wm2="$D_WM2/placement.tsv"
