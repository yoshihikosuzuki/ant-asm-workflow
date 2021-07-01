#!/bin/bash
#SBATCH -J busco
#SBATCH -o busco.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 128
#SBATCH --mem=500G
#SBATCH -t 24:00:00
shopt -s expand_aliases && source ~/.bashrc && set -e || exit 1

ASM=
N_THREADS=128

_ASM=$(basename ${ASM} .gz)
OUT_PREFIX=${_ASM%.*}
OUT_DIR=${OUT_PREFIX}.busco

ml BUSCO
## Metaeuk
busco --download_path $HOME/busco_downloads -c ${N_THREADS} -m genome -l hymenoptera_odb10 -i ${ASM} -o ${OUT_DIR}
## Augustus
busco --download_path $HOME/busco_downloads -c ${N_THREADS} -m genome -l hymenoptera_odb10 --augustus -i ${ASM} -o ${OUT_DIR}.augustus
