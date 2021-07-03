#!/bin/bash
#SBATCH -J asset
#SBATCH -o asset.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=500G
#SBATCH -t 24:00:00
shopt -s expand_aliases && source ~/.bashrc && set -e || exit 1

SCAF=contigs.fasta
READS=hifi.fastq
PB_BAM=contigs.hifi.winnowmap.sorted.bam

OUT_PREFIX=${SCAF%.*}
BED_SCAF=${OUT_PREFIX}.bed
BED_SCAF_GAP=${OUT_PREFIX}.gaps.bed
BED_SCAF_ACC=${OUT_PREFIX}.acc.bed
BED_SCAF_OK=${OUT_PREFIX}.reliable.bed
BED_SCAF_NG=${OUT_PREFIX}.unreliable.bed
PB_SAM=${PB_BAM/.bam/.sam}
PB_PAF=${PB_BAM/.bam/.paf}
BED_PB=${OUT_PREFIX}.pb.bed
BED_PB_OK=${OUT_PREFIX}.pb.reliable.bed

ml asset bedtools samtools

# Contiguous regions and gaps
#samtools faix ${SCAF}
awk '{print $1"\t0\t"$2}' ${SCAF}.fai >${BED_SCAF}
detgaps ${SCAF} >${BED_SCAF_GAP}

bed_to_support() {
    IN_BED=$1
    PREFIX=${IN_BED%.bed}
    BED_NG=${PREFIX}.unreliable.bed
    BED_OK=${PREFIX}.reliable.bed
    bedtools subtract -a ${BED_SCAF} -b ${IN_BED} |
        bedtools merge -d 100 -i - >${BED_NG}
    bedtools subtract -a ${BED_SCAF} -b ${BED_NG} >${BED_OK}
}

# PacBio
samtools view -h ${PB_BAM} -o ${PB_SAM}
paftools.js sam2paf ${PB_SAM} >${PB_PAF}
ast_pb -m3 ${PB_PAF} >${BED_PB}
bed_to_support ${BED_PB}

# Reliable block
acc ${BED_SCAF_GAP} ${BED_PB_OK} |
    awk '$4>=1' |
    bedtools merge -i - >${BED_SCAF_ACC}
bedtools subtract -a ${BED_SCAF} -b ${BED_SCAF_ACC} |
    bedtools merge -d 100 -i - >${BED_SCAF_NG}
bedtools subtract -a ${BED_SCAF} -b ${BED_SCAF_NG} >${BED_SCAF_OK}
info echo -n "Reliable block N50 length = "
awk '{print $3 - $2}' ${BED_SCAF_OK} |
    sort -nr |
    awk '{ sum += $0; print $0, sum }' |
    tac |
    awk 'NR==1 { halftot=$2/2 } lastsize>halftot && $2<halftot { print $1 } { lastsize=$2 }'
