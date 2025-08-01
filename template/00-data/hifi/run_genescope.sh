#!/bin/bash
#SBATCH -J genescope
#SBATCH -o genescope.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=100G
#SBATCH -t 24:00:00
source ../../config.sh
set -eu
ml ${_GENESCOPE} ${_MERQURYFK}
set -x

FASTK_PREFIX=hifi.fastk
K=${HIFI_K}
PLOIDY=${PLOIDY}
HIST_MAX=1000
THRES_ERROR=${HIFI_THRES_ERROR}
N_THREAD=16

OUT_PREFIX=${FASTK_PREFIX/.fastk/.genescope}

Histex -G ${FASTK_PREFIX} -h${HIST_MAX} |
    GeneScopeFK.R -o ${OUT_PREFIX} -p${PLOIDY} -k ${K}
KatGC -pdf -T${N_THREAD} -o${FASTK_PREFIX}.katgc ${FASTK_PREFIX}
PloidyPlot -v -pdf -T${N_THREAD} -e${THRES_ERROR} -o${FASTK_PREFIX}.ploidyplot ${FASTK_PREFIX}
