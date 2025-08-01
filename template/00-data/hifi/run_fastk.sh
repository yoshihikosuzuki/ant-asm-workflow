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

IN_FNAME=hifi.fastq${HIFI_GZ}
K=${HIFI_K}
N_THREAD=16
N_MEMORY=32

_IN_FNAME=$(basename ${IN_FNAME} .gz | sed 's/\.[^.]*$//')
OUT_PREFIX=${_IN_FNAME}.fastk

FastK -k${K} -T${N_THREAD} -M${N_MEMORY} -v -t1 -p -P${TMPDIR} -N${OUT_PREFIX} ${IN_FNAME}
Tabex ${OUT_PREFIX} CHECK
