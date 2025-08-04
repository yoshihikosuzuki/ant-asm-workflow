#!/bin/bash
#SBATCH -J busco
#SBATCH -o busco.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 128
#SBATCH --mem=500G
#SBATCH -t 48:00:00
source ../../../config.sh
set -eu
module load ${_BUSCO}
set -x

ASM=contigs.fasta
BUSCO_DB=${BUSCO_DB}
USE_AUGUSTUS=${USE_AUGUSTUS}
N_CORE=128

OUT_PREFIX=$(basename ${ASM} .gz)
OUT_DIR=${OUT_PREFIX%.*}.busco
if [ "${USE_AUGUSTUS}" = "true" ]; then
    AUGUSTUS_OPTION="--augustus"
else
    AUGUSTUS_OPTION=""
fi

busco -f --update-data -c ${N_CORE} -m genome -l ${BUSCO_DB} -i ${ASM} -o ${OUT_DIR} ${AUGUSTUS_OPTION}

if [ "$AUTO_DEL" = "true" ]; then
    source ./remove_tmp_files.sh
fi
