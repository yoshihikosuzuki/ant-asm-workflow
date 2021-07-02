#!/bin/bash
#SBATCH -J make_index
#SBATCH -o make_index.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=500G
#SBATCH -t 24:00:00
shopt -s expand_aliases && source ~/.bashrc && set -e || exit 1

IN_REF=contigs.fasta

ml samtools
samtools faidx ${IN_REF}

ml bwa
bwa index -p ${IN_REF} ${IN_REF}

#ml blast+
#makeblastdb -in ${IN_REF} -out ${IN_REF} -dbtype nucl -parse_seqids

#ml LAST
#lastdb -uNEAR -R01 ${IN_REF} ${IN_REF}
