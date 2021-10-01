#!/bin/bash
#SBATCH -J fastk
#SBATCH -o fastk.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=40G
#SBATCH -t 24:00:00
shopt -s expand_aliases && source ~/.bashrc && set -e || exit 1
source ../../config.sh

IN_FNAMES="omnic_R1_001.fastq omnic_R2_001.fastq"
K=${OMNIC_K}
N_THREAD=16
N_MEMORY=16
TMP_DIR=tmp

OUT_PREFIX=omnic.fastk

ml Other/FASTK

mkdir -p ${TMP_DIR}
FastK -k${K} -T${N_THREAD} -M${N_MEMORY} -v -t1 -p -P${TMP_DIR} -N${OUT_PREFIX} ${IN_FNAMES}
Tabex ${OUT_PREFIX} CHECK
