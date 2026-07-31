#!/usr/bin/env bash
source /bio/lmod/lmod/init/bash
module purge
module use /bio/package/.modulefiles
module use /nfs/data05/yoshihiko_s/app/.modulefiles

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WR="$REPO/work_real_ont"
SCRIPTS="$WR/scripts"
DAZZ_DB_ROOT="$REPO/DAZZ_DB"
FASTK_ROOT="$REPO/FASTK"
export PATH="$DAZZ_DB_ROOT:$FASTK_ROOT:$REPO:$PATH"

_SEQKIT="seqkit/2.6.0"
_DALIGNER="daligner/2024.01.19"
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

D_REF="$WR/ref"
D_READS="$WR/reads"
D_RAW="$D_READS/raw"
D_TRIO="$WR/trio"
D_EVAL="$WR/eval"
D_HIMAPPER="$WR/himapper"
D_MM2="$WR/mm2"
D_WM="$WR/wm"
D_WM2="$WR/wm2"
mkdir -p "$D_REF" "$D_READS" "$D_RAW" "$D_TRIO" "$D_EVAL" \
  "$D_HIMAPPER" "$D_MM2" "$D_WM" "$D_WM2"

# HG002 R10.4.1 ultra-long reads, 2022-11-16 runs, Guppy 6.3.8 SUP.
ONT_FASTQ_BASE="https://human-pangenomics.s3.us-west-2.amazonaws.com/T2T/scratch/HG002/sequencing/ont/11_16_22_R1041_HG002_UL_Kit14_260_Guppy_6.3.8"
ONT_FASTQ_NAMES=(
  "11_16_22_R1041_HG002_UL_Kit14_260_4_Guppy_6.3.8_sup.pass.lt100kb.fastq.gz"
  "11_16_22_R1041_HG002_UL_Kit14_260_5_Guppy_6.3.8_sup.pass.lt100kb.fastq.gz"
  "11_16_22_R1041_HG002_UL_Kit14_260_6_Guppy_6.3.8_sup.pass.lt100kb.fastq.gz"
)
ONT_FASTQ_SIZES=(32818692274 33611072332 31410844872)
ONT_FASTQS=(
  "$D_RAW/${ONT_FASTQ_NAMES[0]}"
  "$D_RAW/${ONT_FASTQ_NAMES[1]}"
  "$D_RAW/${ONT_FASTQ_NAMES[2]}"
)

# Later candidate: HG002 R10.4.1 ultra-long reads, 2024-09-16,
# Dorado 0.9.1, dna_r10.4.1_e8.2_400bps_sup@v5.0.0.
# ONT_BAM_BASE="https://human-pangenomics.s3.us-west-2.amazonaws.com/submissions/5b73fa0e-658a-4248-b2b8-cd16155bc157--UCSC_GIAB_R1041_nanopore/HG002_R1041_UL/dorado/supv5.0.0_wMods"
# ONT_BAM_NAMES=(
#   "09_16_24_R1041_UL_HG002_7C_dorado0.9.1_sup5.0.0_5mCG_5hmCG.bam"
#   "09_16_24_R1041_UL_HG002_7T_dorado0.9.1_sup5.0.0_5mCG_5hmCG.bam"
# )
# ONT_BAM_SIZES=(57463442321 52349634898)
# ONT_BAMS=(
#   "$D_RAW/${ONT_BAM_NAMES[0]}"
#   "$D_RAW/${ONT_BAM_NAMES[1]}"
# )

REF_ALL="$D_REF/dip.fa"
REF="$D_REF/chr6.fa"
DIP_REPETITIVE="$D_REF/dip.repetitive.txt"
CHR6_MAT="chr6_MATERNAL"
CHR6_PAT="chr6_PATERNAL"

MIN_READ_LENGTH=50000
MIN_QUERY_COVERAGE=0.90
MIN_TRUTH_MAPQ=60

TRUTH_SOURCE_PAF="$D_READS/truth_source.wm2.paf"
READS="$D_READS/chr6.ont-ul.guppy6.3.8-sup.gt50kb.reads.fa"
READS_BASE="$(basename "$READS" .fa)"
READ_NAMES="$D_READS/$READS_BASE.names"
READS_MAPPING_TRUTH="$D_READS/$READS_BASE.mapping_truth.tsv"
READS_CONFIDENT_MAPPING_TRUTH="$D_READS/$READS_BASE.confident.mapping_truth.tsv"

# ONT mainline defaults.  Keep seed density at w=s=100; e and the unique-seed
# selection weight were selected with the chain-aware chr6 ONT comparison.
KMER=32
HIMAPPER_WS=100
HIMAPPER_E=0.70
HIMAPPER_TRACE_SPACING=100
DB_BLOCK_SIZE=200

CN_UNIQ_W=200
CN_TIER_SRC=1
CN_HASH_TIE=1
CN_DUMP=0

ONT_TMP_ROOT=${ONT_TMP_ROOT:-/grid2/yoshihiko_s/himapper_ont}
