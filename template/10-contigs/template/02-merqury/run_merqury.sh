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

READS=hifi.fastq
CONTIGS=contigs.fasta

# NOTE: The value of `K` should be determined by the following command:
#GENOME_SIZE=300000000
#best_k.sh ${GENOME_SIZE}
K=19

N_THREADS=128
N_MEMORY=500

READS_MERYL=${READS}.meryl
_READS=$(basename ${READS} .gz)
_CONTIGS=$(basename ${CONTIGS} .gz)
OUT_PREFIX=${_CONTIGS%.*}.${_READS%.*}.merqury

ml Other/merqury

meryl count k=${K} memory=${N_MEMORY} threads=${N_THREADS} output ${READS_MERYL} ${READS}
merqury.sh ${READS_MERYL} ${CONTIGS} ${OUT_PREFIX}
