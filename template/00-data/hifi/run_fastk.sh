#!/bin/bash
#SBATCH -J fastk
#SBATCH -o fastk.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=500G
#SBATCH -t 24:00:00
shopt -s expand_aliases && source ~/.bashrc && set -e || exit 1

IN_FNAME=hifi.fastq
K=40
N_THREAD=16
N_MEMORY=16
TMP_DIR=tmp

OUT_PREFIX=${IN_FNAME%.*}.fastk

ml Other/FASTK
mkdir -p ${TMP_DIR}

FastK -k${K} -T${N_THREAD} -M${N_MEMORY} -v -t1 -p -P${TMP_DIR} -N${OUT_PREFIX} ${IN_FNAME}
Tabex ${OUT_PREFIX} CHECK
