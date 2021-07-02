#!/bin/bash
#SBATCH -J winnowmap
#SBATCH -o winnowmap.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 128
#SBATCH --mem=500G
#SBATCH -t 24:00:00
shopt -s expand_aliases && source ~/.bashrc && set -e || exit 1

# NOTE: This script is specifically for purge_dups, which generate reads.paf and contigs.paf

CONTIGS=contigs.fasta
READS=hifi.fastq
K=15
N_THREADS=128

# NOTE: Do not overload `meryl` after this
ml samtools winnowmap purge_dups

## Contig vs read

REF=${CONTIGS}
_REF=$(basename ${REF} .gz)
_READS=$(basename ${READS} .gz)
OUT_MERYL=${_REF}.winnowmap_meryl_k${K}
OUT_REP=${_REF}.winnowmap_rep_k${K}

meryl count k=${K} output ${OUT_MERYL} ${REF}
meryl print greater-than distinct=0.9998 ${OUT_MERYL} >${OUT_REP}

winnowmap -t${N_THREADS} -xmap-pb -W ${OUT_REP} ${REF} ${READS} >reads.paf

## Contig vs contig

CONTIGS_PREFIX=${CONTIGS%.gz}
CONTIGS_PREFIX=${CONTIGS_PREFIX%.*}
CONTIGS_SPLIT=${CONTIGS_PREFIX}.split.fa

split_fa ${CONTIGS} >${CONTIGS_SPLIT}

REF=${CONTIGS_SPLIT}
QUERY=${CONTIGS_SPLIT}
_REF=$(basename ${REF} .gz)
_QUERY=$(basename ${QUERY} .gz)
OUT_MERYL=${_REF}.winnowmap_meryl_k${K}
OUT_REP=${_REF}.winnowmap_rep_k${K}

meryl count k=${K} output ${OUT_MERYL} ${REF}
meryl print greater-than distinct=0.9998 ${OUT_MERYL} >${OUT_REP}
winnowmap -t${N_THREADS} -xasm5 -DP -W ${OUT_REP} ${REF} ${QUERY} >contigs.paf
