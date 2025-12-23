#!/bin/bash
#SBATCH -J compleasm
#SBATCH -o compleasm.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 128
#SBATCH --mem=500G
#SBATCH -t 48:00:00
source ../../../config.sh
set -eu
module load ${_COMPLEASM}
set -x

ASM=scaffolds.fasta
N_CORE=128

OUT_PREFIX=$(basename ${ASM} .gz)
OUT_DIR=${OUT_PREFIX%.*}.compleasm

# compleasm download ${BUSCO_DB}
compleasm run -t ${N_CORE} -l ${BUSCO_DB} -L ${BUSCO_DB_DIR} -a ${ASM} -o ${OUT_DIR}

if [ "$AUTO_DEL" = "true" ]; then
    source ./remove_tmp_files.sh
fi
