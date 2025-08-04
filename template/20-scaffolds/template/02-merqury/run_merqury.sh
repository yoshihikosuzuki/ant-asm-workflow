#!/bin/bash
#SBATCH -J merqury
#SBATCH -o merqury.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 128
#SBATCH --mem=500G
#SBATCH -t 24:00:00
source ../../../config.sh
set -eu
module load ${_MERQURY}
set -x

SCAFS=scaffolds.fasta
READS=hifi.fastq${HIFI_GZ}
K=${MERQURY_K}
N_THREADS=128
N_MEMORY=500

READS_MERYL=${READS}.meryl
_READS=$(basename ${READS} .gz | sed 's/\.[^.]*$//')
_SCAFS=$(basename ${SCAFS} .gz | sed 's/\.[^.]*$//')
OUT_PREFIX=${_SCAFS}.${_READS}.merqury

meryl count k=${K} memory=${N_MEMORY} threads=${N_THREADS} output ${READS_MERYL} ${READS}
merqury.sh ${READS_MERYL} ${SCAFS} ${OUT_PREFIX}

if [ "$AUTO_DEL" = "true" ]; then
    source ./remove_tmp_files.sh
fi
