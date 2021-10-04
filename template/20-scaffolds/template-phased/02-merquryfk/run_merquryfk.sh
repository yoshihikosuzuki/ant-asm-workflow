#!/bin/bash
#SBATCH -J merquryfk
#SBATCH -o merquryfk.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=40G
#SBATCH -t 24:00:00
shopt -s expand_aliases && source ~/.bashrc && set -e || exit 1
source ../../../config.sh

READS=hifi.fastq
HAP1=scaffolds.hap1.fasta
HAP2=scaffolds.hap2.fasta
K=${MERQURY_K}
N_THREADS=16
TMP_DIR=tmp

_READS=$(basename ${READS} .gz)
_HAP1=$(basename ${HAP1} .gz)
_HAP2=$(basename ${HAP2} .gz)
_READS=${_READS%.*}
_HAP1=${_HAP1%.*}
_HAP2=${_HAP2%.*}
READS_FASTK=${_READS}.fastk
HAP1_FASTK=${_HAP1}.fastk
HAP2_FASTK=${_HAP2}.fastk
OUT_PREFIX=${_HAP1#scaffolds.}.${_HAP2#scaffolds.}.${_READS}.merqury

ml ${_MERQURYFK}

mkdir -p ${TMP_DIR}
FastK -k${K} -T${N_THREADS} -v -t1 -P${TMP_DIR} -N${READS_FASTK} ${READS}
FastK -k${K} -T${N_THREADS} -v -t1 -p -P${TMP_DIR} -N${HAP1_FASTK} ${HAP1}
FastK -k${K} -T${N_THREADS} -v -p:${READS_FASTK} -N${HAP1_FASTK}.${READS_FASTK} ${HAP1}
FastK -k${K} -T${N_THREADS} -v -t1 -p -P${TMP_DIR} -N${HAP2_FASTK} ${HAP2}
FastK -k${K} -T${N_THREADS} -v -p:${READS_FASTK} -N${HAP2_FASTK}.${READS_FASTK} ${HAP2}

CNspectra -v -pdf -T${N_THREADS} ${READS_FASTK} ${HAP1_FASTK} ${HAP2_FASTK} ${OUT_PREFIX}

if [ "$AUTO_DEL" = "true" ]; then
    source ./remove_tmp_files.sh
fi
