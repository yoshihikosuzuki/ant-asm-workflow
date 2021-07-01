#!/bin/bash
#SBATCH -J bwa
#SBATCH -o bwa.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 128
#SBATCH --mem=500G
#SBATCH -t 24:00:00

REF_FASTA=
QUERY_FASTA=
OUT_BAM=

N_THREADS=128

ml samtools bwa

bwa index ${IN_REF}
bwa mem -t${N_THREADS} ${IN_REF} ${QUERY_FASTA} |
    samtools sort -@${N_THREADS} -o ${OUT_BAM}
samtools index -@${N_THREADS} ${OUT_BAM}
