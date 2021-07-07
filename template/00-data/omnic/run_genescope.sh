#!/bin/bash
#SBATCH -J genescope
#SBATCH -o genescope.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=500G
#SBATCH -t 24:00:00
shopt -s expand_aliases && source ~/.bashrc && set -e || exit 1

FASTK_PREFIX=omnic.fastk
K=21
PLOIDY=2
HIST_MAX=1000

OUT_PREFIX=${FASTK_PREFIX%.fastk}.genescope

ml Other/genescope

Histex -G ${FASTK_PREFIX} -h${HIST_MAX} | GeneScopeFK.R -o ${OUT_PREFIX} -p${PLOIDY} -k ${K}
