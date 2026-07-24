#!/usr/bin/env bash
source /bio/lmod/lmod/init/bash
module purge
module use /bio/package/.modulefiles
module use /nfs/data05/yoshihiko_s/app/.modulefiles
export TMPDIR=/grid2/yoshihiko_s/
mkdir -p "$TMPDIR"

# Repository
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WR="$REPO/work_real"

# Tools
SCRIPTS="$WR/scripts"
export PATH="$REPO:$PATH"

_SEQKIT="seqkit/2.6.0"
_DAZZ_DB="dazz_db/2025.02.17"
_DALIGNER="daligner/2024.01.19"
_FASTK="fastk/2025.09.13"
_SAMTOOLS="samtools/1.18"
_MINIMAP2="minimap2/2.26"
_WINNOWMAP="winnowmap/2.03"
_MERYL="meryl/1.3"
_YAK="yak/0.1"
_PYTHON="python/3.10.14"

TIME_BIN="/usr/bin/time"
run_timed() {
  local time_file=$1
  shift
  "$TIME_BIN" -v -o "$time_file" "$@"
}

plot_coverage() {
  local output=$1
  shift
  python3 "$SCRIPTS/plot_coverage.py" "$output" "$@"
}

# Work directories
D_REF="$WR/ref"
D_READS="$WR/reads"
D_TRIO="$WR/trio"
D_REGIONS="$WR/regions"
D_EVAL="$WR/eval"
D_HIMAPPER="$WR/himapper"
D_MM2="$WR/mm2"
D_WM2="$WR/wm2"

mkdir -p "$D_REF" "$D_READS" "$D_TRIO" "$D_REGIONS" "$D_EVAL" \
  "$D_HIMAPPER" "$D_MM2" "$D_WM2"

# Whole-genome inputs
REF_ALL="$D_REF/dip.fa"
READS_ALL="$D_READS/hg002.hifi.fastq.gz"

# chr6 datasets
CHR6_MAT="chr6_MATERNAL"
CHR6_PAT="chr6_PATERNAL"
REF="$D_REF/chr6.fa"
READS="$D_READS/chr6.reads.fa"
READS_BASE="$(basename "$READS" .fa)"
READ_NAMES="$D_READS/$READS_BASE.names"
READS_MAPPING_TRUTH="$D_READS/$READS_BASE.mapping_truth.tsv"
READS_CONFIDENT_MAPPING_TRUTH="$D_READS/$READS_BASE.confident.mapping_truth.tsv"

# Parameters
KMER=32
