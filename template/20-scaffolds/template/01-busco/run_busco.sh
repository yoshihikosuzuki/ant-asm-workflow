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

ASM=scaffolds.fasta
LINEAGE=$BUSCO_DB
DOWNLOAD_PATH="$HOME/busco_downloads"
N_THREADS=16

OUT_PREFIX=$(basename ${ASM} .gz)
OUT_DIR=${OUT_PREFIX%.*}.busco
AUGUSTUS_PATH=${OUT_DIR}_augustus

ml ${_BUSCO}

## Case 1. Using Metaeuk for gene annotation
#busco -f --download_path ${DOWNLOAD_PATH} -c ${N_THREADS} -m genome -l ${LINEAGE} -i ${ASM} -o ${OUT_DIR}

## Case 2. Using Augustus for gene annotation
busco -f --download_path ${DOWNLOAD_PATH} -c ${N_THREADS} -m genome -l ${LINEAGE} --augustus -i ${ASM} -o ${AUGUSTUS_PATH}

if [ "$AUTO_DEL" = "true" ]; then
    source ./remove_tmp_files.sh
fi
