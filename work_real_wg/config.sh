#!/usr/bin/env bash
source /bio/lmod/lmod/init/bash
module purge
module use /bio/package/.modulefiles
module use /nfs/data05/yoshihiko_s/app/.modulefiles

# Repository
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WR="$REPO/work_real_wg"
SCRIPTS="$WR/scripts"
DAZZ_DB_ROOT="$REPO/DAZZ_DB"
FASTK_ROOT="$REPO/FASTK"
export PATH="$DAZZ_DB_ROOT:$FASTK_ROOT:$REPO:$PATH"

# Tools
_SEQKIT="seqkit/2.6.0"
_DALIGNER="daligner/2024.01.19"
_MINIMAP2="minimap2/2.26"
_YAK="yak/0.1"
_PYTHON="python/3.10.14"

TIME_BIN="/usr/bin/time"
run_timed() {
  local time_file=$1
  shift
  "$TIME_BIN" -v -o "$time_file" "$@"
}

# Work directories
D_REF="$WR/ref"
D_READS="$WR/reads"
D_TRIO="$WR/trio"
D_EVAL="$WR/eval"
D_HIMAPPER="$WR/himapper"
D_MM2="$WR/mm2"
mkdir -p "$D_REF" "$D_READS" "$D_TRIO" "$D_EVAL" "$D_HIMAPPER" "$D_MM2"

# Evaluation annotations and derived exclusion mask.
EVAL_ANNOTATION_DIR="$WR/experiments/hotspot_coverage/annotations"
EVAL_RDNA_GAP_BED="$EVAL_ANNOTATION_DIR/hg002v1.1.rDNA_GAP.envelope.bed"
EVAL_CENSAT_BED="$EVAL_ANNOTATION_DIR/hg002v1.1.cenSatv2.0.bed"
EVAL_EXCLUSION_BED="$D_EVAL/hg002v1.1.rDNA_GAP_and_active_HOR.bed"
EVAL_EXCLUSION_METADATA="$D_EVAL/hg002v1.1.rDNA_GAP_and_active_HOR.metadata.tsv"
EVAL_ALIGNMENT_STRATA="$D_EVAL/alignment_strata.tsv"

# Inputs and derived read cohort
REF="$D_REF/dip.fa"
READS_ALL="$D_READS/hg002.hifi.fastq.gz"
TRUTH_SOURCE_PAF="$D_READS/truth_source.wm2.paf"
READS="$D_READS/wg.reads.fa"
READS_BASE="wg.reads"
READ_NAMES="$D_READS/$READS_BASE.names"
READS_MAPPING_TRUTH="$D_READS/$READS_BASE.mapping_truth.tsv"
READS_CONFIDENT_MAPPING_TRUTH="$D_READS/$READS_BASE.confident.mapping_truth.tsv"

# Fixed cohort and mapper parameters
MIN_READ_LENGTH=10000
MIN_QUERY_COVERAGE=0.90
MIN_TRUTH_MAPQ=60
KMER=32
HIMAPPER_WS=${HIMAPPER_WS:-100}
MM2_K=19
MM2_W=19
DB_BLOCK_SIZE=200

# himapper candidate selection
CN_UNIQ_W=11.5
CN_TIER_SRC=1
CN_HASH_TIE=1
CN_DUMP=0

# Node-local/intermediate scratch; override at submission if needed.
WG_TMP_ROOT=${WG_TMP_ROOT:-/grid2/yoshihiko_s/himapper_wg}
