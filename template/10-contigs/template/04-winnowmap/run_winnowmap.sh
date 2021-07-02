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

REF=contigs.fasta
READS=hifi.fastq
N_THREADS=128

# NOTE: Do not overload `meryl` after this
ml samtools winnowmap

TYPE="map-pb"
K=15

_REF=$(basename ${REF} .gz)
_READS=$(basename ${READS} .gz)
OUT_MERYL=${_REF}.winnowmap_meryl_k${K}
OUT_REP=${_REF}.winnowmap_rep_k${K}
OUT_PREFIX=${_REF%.*}.${_READS%.*}.winnowmap
OUT_BAM=${OUT_PREFIX}.sorted.bam

meryl count k=${K} output ${OUT_MERYL} ${REF}
meryl print greater-than distinct=0.9998 ${OUT_MERYL} >${OUT_REP}
## With secondary mappings
#winnowmap -t${N_THREADS} -ax ${TYPE} -W ${OUT_REP} ${REF} ${READS} |
#    samtools sort -@${N_THREADS} -o ${OUT_BAM}
## Without secondary mappings
winnowmap -t${N_THREADS} -ax ${TYPE} --secondary=no --eqx -Y -W ${OUT_REP} ${REF} ${READS} |
    samtools sort -@${N_THREADS} -o ${OUT_BAM}
samtools index -@${N_THREADS} ${OUT_BAM}
