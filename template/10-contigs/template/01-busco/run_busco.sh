#!/bin/bash
#SBATCH -J busco
#SBATCH -o busco.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 128
#SBATCH --mem=500G
#SBATCH -t 48:00:00
shopt -s expand_aliases && source ~/.bashrc && set -e || exit 1

ASM=contigs.fasta
LINEAGE="hymenoptera_odb10"
DOWNLOAD_PATH="$HOME/busco_downloads"
N_THREADS=128

OUT_PREFIX=$(basename ${ASM} .gz)
OUT_PREFIX=${OUT_PREFIX%.*}
OUT_DIR=${OUT_PREFIX}.busco

ml BUSCO

## Case 1. Using Metaeuk for gene annotation
#busco --download_path ${DOWNLOAD_PATH} -c ${N_THREADS} -m genome -l ${LINEAGE} -i ${ASM} -o ${OUT_DIR}

## Case 2. Using Augustus for gene annotation
busco --download_path ${DOWNLOAD_PATH} -c ${N_THREADS} -m genome -l ${LINEAGE} --augustus -i ${ASM} -o ${OUT_DIR}_augustus
