#!/bin/bash
#SBATCH -J merquryfk
#SBATCH -o merquryfk.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=40G
#SBATCH -t 24:00:00
source ../../../config.sh
set -eu
module load ${_MERQURYFK}
set -x

SCAFS=scaffolds.fasta
READS=hifi.fastq${HIFI_GZ}
K=${MERQURY_K}
N_THREADS=16

_READS=$(basename ${READS} .gz | sed 's/\.[^.]*$//')
_SCAFS=$(basename ${SCAFS} .gz | sed 's/\.[^.]*$//')
READS_FASTK=${_READS}.fastk
SCAFS_FASTK=${_SCAFS}.fastk
OUT_PREFIX=${_SCAFS}.${_READS}.merqury

FastK -k${K} -T${N_THREADS} -v -t1 -P${TMPDIR} -N${READS_FASTK} ${READS}
FastK -k${K} -T${N_THREADS} -v -t1 -p -P${TMPDIR} -N${SCAFS_FASTK} ${SCAFS}
FastK -k${K} -T${N_THREADS} -v -p:${READS_FASTK} -N${SCAFS_FASTK}.${READS_FASTK} ${SCAFS}

CNspectra -v -pdf -T${N_THREADS} ${READS_FASTK} ${SCAFS_FASTK} ${OUT_PREFIX}

if [ "$AUTO_DEL" = "true" ]; then
    source ./remove_tmp_files.sh
fi
