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

CONTIGS=contigs.fasta
READS=hifi.fastq${HIFI_GZ}
K=${MERQURY_K}
N_THREADS=16

_READS=$(basename ${READS} .gz | sed 's/\.[^.]*$//')
_CONTIGS=$(basename ${CONTIGS} .gz | sed 's/\.[^.]*$//')
READS_FASTK=${_READS}.fastk
CONTIGS_FASTK=${_CONTIGS}.fastk
OUT_PREFIX=${_CONTIGS}.${_READS}.merqury

FastK -k${K} -T${N_THREADS} -v -t1 -P${TMPDIR} -N${READS_FASTK} ${READS}
FastK -k${K} -T${N_THREADS} -v -t1 -p -P${TMPDIR} -N${CONTIGS_FASTK} ${CONTIGS}
FastK -k${K} -T${N_THREADS} -v -p:${READS_FASTK} -N${CONTIGS_FASTK}.${READS_FASTK} ${CONTIGS}

CNspectra -v -pdf -T${N_THREADS} ${READS_FASTK} ${CONTIGS_FASTK} ${OUT_PREFIX}

if [ "$AUTO_DEL" = "true" ]; then
    source ./remove_tmp_files.sh
fi
