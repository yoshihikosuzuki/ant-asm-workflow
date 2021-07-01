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

N_THREADS=128
N_MEMORY=500

ml merqury

#GENOME_SIZE=
#best_k.sh ${GENOME_SIZE}
K=

## Trio assembly

CHILD_READS=
MAT_READS=
PAT_READS=
MAT_CONTIGS=
PAT_CONTIGS=
OUT_PREFIX=

meryl count k=${K} memory=${N_MEMORY} threads=${N_THREADS} output ${CHILD_READS}.meryl ${CHILD_READS}
meryl count k=${K} memory=${N_MEMORY} threads=${N_THREADS} output ${MAT_READS}.meryl ${MAT_READS}
meryl count k=${K} memory=${N_MEMORY} threads=${N_THREADS} output ${PAT_READS}.meryl ${PAT_READS}
$MERQURY/trio/hapmers.sh ${MAT_READS}.meryl ${PAT_READS}.meryl ${CHILD_READS}.meryl
merqury.sh ${CHILD_READS}.meryl ${MAT_READS}.hapmer.meryl ${PAT_READS}.hapmer.meryl ${MAT_CONTIGS} ${PAT_CONTIGS} ${OUT_PREFIX}

## Merged assembly

READS=
CONTIGS=
OUT_PREFIX=

meryl count k=${K} memory=${N_MEMORY} threads=${N_THREADS} output ${READS}.meryl ${READS}
merqury.sh ${READS}.meryl ${CONTIGS} ${OUT_PREFIX}

## Haplotype-resolved assembly

READS=
PRIMARY_CONTIGS=
ALTERNATE_CONTIGS=
OUT_PREFIX=

meryl count k=${K} memory=${N_MEMORY} threads=${N_THREADS} output ${READS}.meryl ${READS}
merqury.sh ${READS}.meryl ${PRIMARY_CONTIGS} ${ALTERNATE_CONTIGS} ${OUT_PREFIX}
