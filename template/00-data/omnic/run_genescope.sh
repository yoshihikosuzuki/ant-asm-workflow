#!/bin/bash
#SBATCH -J genescope
#SBATCH -o genescope.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=100G
#SBATCH -t 24:00:00
shopt -s expand_aliases && source ~/.bashrc && set -e || exit 1
source ../../config.sh

FASTK_PREFIX=omnic.fastk
K=${OMNIC_K}
PLOIDY=${PLOIDY}
HIST_MAX=1000
THRES_ERROR=${OMNIC_THRES_ERROR}
N_THREAD=16

OUT_PREFIX=${FASTK_PREFIX/.fastk/.genescope}

ml Other/genescope Other/MerquryFK

Histex -G ${FASTK_PREFIX} -h${HIST_MAX} |
    GeneScopeFK.R -o ${OUT_PREFIX} -p${PLOIDY} -k ${K}
KatGC -pdf -T${N_THREAD} -o${FASTK_PREFIX}.katgc ${FASTK_PREFIX}
PloidyPlot -v -pdf -T${N_THREAD} -e${ERROR_THRES} -o${FASTK_PREFIX}.ploidyplot ${FASTK_PREFIX}
