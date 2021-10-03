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
SCAFS=scaffolds.fasta
K=${MERQURY_K}
N_THREADS=16
TMP_DIR=tmp

_READS=$(basename ${READS} .gz)
_SCAFS=$(basename ${SCAFS} .gz)
_READS=${_READS%.*}
_SCAFS=${_SCAFSS%.*}
READS_FASTK=${_READS}.fastk
SCAFS_FASTK=${_SCAFS}.fastk
OUT_PREFIX=${_SCAFS}.${_READS}.merqury

ml ${_MERQURYFK}

mkdir -p ${TMP_DIR}
FastK -k${K} -T${N_THREADS} -v -t1 -P${TMP_DIR} -N${READS_FASTK} ${READS}
FastK -k${K} -T${N_THREADS} -v -t1 -p -P${TMP_DIR} -N${SCAFS_FASTK} ${SCAFS}
FastK -k${K} -T${N_THREADS} -v -p:${READS_FASTK} -N${SCAFS_FASTK}.${READS_FASTK} ${SCAFS}

CNspectra -v -pdf -T${N_THREADS} ${READS_FASTK} ${SCAFS_FASTK} ${OUT_PREFIX}

if [ "$AUTO_DEL" = "true" ]; then
    source ./remove_tmp_files.sh
fi
