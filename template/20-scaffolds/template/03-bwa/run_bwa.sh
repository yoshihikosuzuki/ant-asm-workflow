#!/bin/bash
#SBATCH -J bwa
#SBATCH -o bwa.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 128
#SBATCH --mem=500G
#SBATCH -t 24:00:00
shopt -s expand_aliases && source ~/.bashrc && set -e || exit 1

SCAF=scaffolds.fasta
HIC_READS_1=omnic_R1_001.fastq
HIC_READS_2=omnic_R2_001.fastq
N_THREADS=128

_SCAF=$(basename ${SCAF} .gz)
_READS=$(basename ${HIC_READS_1} .gz)
_READS=${_READS%_R*}
OUT_BAM=${_SCAF%.*}.${_READS%.*}.sorted.bam

ml samtools bwa

#bwa index ${SCAF}
bwa mem -t${N_THREADS} -5SP -B8 ${SCAF} ${HIC_READS_1} ${HIC_READS_2} |
    samtools sort -@${N_THREADS} -o ${OUT_BAM}
samtools index -@${N_THREADS} ${OUT_BAM}

# Coverage
BIN_SIZE=1000
N_THREADS=4
OUT_PREFIX=${OUT_BAM}

ml Other/mosdepth

mosdepth -t${N_THREADS} -b ${BIN_SIZE} -n -x ${OUT_PREFIX} ${OUT_BAM}
zcat ${OUT_PREFIX}.regions.bed.gz > ${OUT_PREFIX}.regions.bedgraph
