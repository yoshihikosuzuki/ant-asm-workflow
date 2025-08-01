#!/bin/bash
#SBATCH -J fastk
#SBATCH -o fastk.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=40G
#SBATCH -t 24:00:00
source ../../config.sh
set -eu
ml ${_FASTK}
set -x

IN_FNAMES="omnic_R1_001.fastq${OMNIC_GZ} omnic_R2_001.fastq${OMNIC_GZ}"
K=${OMNIC_K}
N_THREAD=16
N_MEMORY=32

OUT_PREFIX=omnic.fastk

FastK -k${K} -T${N_THREAD} -M${N_MEMORY} -v -t1 -p -P${TMPDIR} -N${OUT_PREFIX} ${IN_FNAMES}
Tabex ${OUT_PREFIX} CHECK
