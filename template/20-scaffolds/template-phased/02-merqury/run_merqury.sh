#!/bin/bash
#SBATCH -J merqury
#SBATCH -o merqury.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 128
#SBATCH --mem=500G
#SBATCH -t 24:00:00
shopt -s expand_aliases && source ~/.bashrc && set -e || exit 1
source ../../../config.sh

READS=hifi.fastq
HAP1=scaffolds.hap1.fasta
HAP2=scaffolds.hap2.fasta
K=${MERQURY_K}
N_THREADS=128
N_MEMORY=500

READS_MERYL=${READS}.meryl
_READS=$(basename ${READS} .gz)
_HAP1=$(basename ${HAP1} .gz)
_HAP2=$(basename ${HAP2} .gz)
_READS=${_READS%.*}
_HAP1=${_HAP1%.*}
_HAP2=${_HAP2%.*}
OUT_PREFIX=${_HAP1#scaffolds.}.${_HAP2#scaffolds.}.${_READS}.merqury

ml ${_MERQURY}

meryl count k=${K} memory=${N_MEMORY} threads=${N_THREADS} output ${READS_MERYL} ${READS}
merqury.sh ${READS_MERYL} ${HAP1} ${HAP2} ${OUT_PREFIX}

if [ "$AUTO_DEL" = "true" ]; then
    source ./remove_tmp_files.sh
fi
