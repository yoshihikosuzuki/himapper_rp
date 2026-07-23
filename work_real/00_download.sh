#!/usr/bin/env bash
#SBATCH -J hg002_download
#SBATCH -o 00_download.log
#SBATCH -p rl8
#SBATCH -c 8
#SBATCH -n 1
#SBATCH -N 1
set -euo pipefail
source config.sh
module load ${_SAMTOOLS}
set -x
THREADS=${SLURM_CPUS_PER_TASK:-8}

GIAB="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab"
HPGP="https://s3-us-west-2.amazonaws.com/human-pangenomics"

# HG002 HiFi reads
HIFI_DIR="$GIAB/data/AshkenazimTrio/HG002_NA24385_son/PacBio_HiFi-Revio_20231031"
HIFI_BAMS="HG002_PacBio-Revio_m84039_230928_213653_s3.hifi_reads.bam HG002_PacBio-Revio_m84039_231005_222902_s1.hifi_reads.bam"

# HG002 diploid assembly
ASM_MAT_URL="$HPGP/T2T/HG002/assemblies/hg002v1.1.mat.fasta.gz"
ASM_PAT_URL="$HPGP/T2T/HG002/assemblies/hg002v1.1.pat.fasta.gz"

# Trio parent reads
HG003_DIR="$GIAB/data/AshkenazimTrio/HG003_NA24149_father/NIST_Illumina_2x250bps/reads"
HG003_GLOB="D2_S1_*.fastq.gz,D2_S2_L001_*.fastq.gz"
HG004_DIR="$GIAB/data/AshkenazimTrio/HG004_NA24143_mother/NIST_Illumina_2x250bps/reads"
HG004_GLOB="D3_S1_*.fastq.gz,D3_S3_L001_R[12]_00[12].fastq.gz"

# GRCh38 regions
CMRG_BED_URL="$GIAB/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/SupplementaryFiles/GRCh38_CMRG_benchmark_gene_coordinates.bed"
STRAT_BASE="$GIAB/release/genome-stratifications/v3.6/GRCh38@all"
SEGDUP_BED_URL="$STRAT_BASE/SegmentalDuplications/GRCh38_segdups.bed.gz"
MHC_BED_URL="$STRAT_BASE/OtherDifficult/GRCh38_MHC.bed.gz"

dl() {
  local url=$1 dst=$2
  wget -q -O "$dst" "$url"
}

dl_glob() {
  local base=$1 patterns=$2 dst=$3
  mkdir -p "$dst"
  wget -q -r -np -nd -A "$patterns" -P "$dst" "$base/"
}

# Combine HiFi BAMs
bams=()
for bam in $HIFI_BAMS; do
  dl "$HIFI_DIR/$bam" "$D_READS/$bam"
  bams+=("$D_READS/$bam")
done
samtools cat "${bams[@]}" |
  samtools fastq -@ "$THREADS" - |
  bgzip -@ "$THREADS" -c > "$READS_ALL"

# Download the diploid assembly
dl "$ASM_MAT_URL"   "$D_REF/hg002.mat.fa.gz"
dl "$ASM_PAT_URL"   "$D_REF/hg002.pat.fa.gz"

# Download trio parent reads
dl_glob "$HG003_DIR" "$HG003_GLOB" "$D_TRIO/hg003"
dl_glob "$HG004_DIR" "$HG004_GLOB" "$D_TRIO/hg004"

# Download GRCh38 regions
dl "$CMRG_BED_URL"    "$D_REGIONS/CMRG.grch38.bed"
dl "$SEGDUP_BED_URL"  "$D_REGIONS/segdups.grch38.bed.gz"
dl "$MHC_BED_URL"     "$D_REGIONS/MHC.grch38.bed.gz"
