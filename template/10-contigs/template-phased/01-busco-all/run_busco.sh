#!/bin/bash
#SBATCH -J busco
#SBATCH -o busco.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=500G
#SBATCH -t 48:00:00
shopt -s expand_aliases && source ~/.bashrc && set -e || exit 1
source ../../../config.sh

ASM=contigs.fasta
N_THREADS=16

OUT_PREFIX=$(basename ${ASM} .gz)
OUT_DIR=${OUT_PREFIX%.*}.busco

ml ${_BUSCO}

SHARED_ARGS="-f --update-data -c ${N_THREADS} -m genome -l ${BUSCO_DB} -i ${ASM}"
## Case 1. Using Metaeuk for gene annotation
#busco ${SHARED_ARGS} -o ${OUT_DIR}
## Case 2. Using Augustus for gene annotation
busco ${SHARED_ARGS} -o ${OUT_DIR}_augustus --augustus

if [ "$AUTO_DEL" = "true" ]; then
    source ./remove_tmp_files.sh
fi
