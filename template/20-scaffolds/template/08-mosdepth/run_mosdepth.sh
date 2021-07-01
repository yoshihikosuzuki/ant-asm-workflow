#!/bin/bash
#SBATCH -J mosdepth
#SBATCH -o mosdepth.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=500G
#SBATCH -t 10:00:00
shopt -s expand_aliases && source ~/.bashrc && set -e || exit 1

IN_BAM=
BIN_SIZE=1000
N_THREADS=4

OUT_PREFIX=${IN_BAM}

ml samtools mosdepth
#samtools index -@${N_THREADS} ${IN_BAM}
mosdepth -t${N_THREADS} -b ${BIN_SIZE} -n -x ${OUT_PREFIX} ${IN_BAM}
zcat ${OUT_PREFIX}.regions.bed.gz > ${OUT_PREFIX}.regions.bedgraph
