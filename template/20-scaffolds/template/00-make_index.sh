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

IN_REF=scaffolds.fasta

ml samtools
samtools faidx ${IN_REF}

ml bwa
bwa index -p ${IN_REF} ${IN_REF}

ml Other/3d-dna
3d-dna-fasta2assembly ${IN_REF} >${IN_REF%.*}.assembly
awk 'NF == 3 {print substr($1,2) "\t" $3}' ${IN_REF%.*}.assembly >${IN_REF%.*}.chrom_sizes

#ml blast+
#makeblastdb -in ${IN_REF} -out ${IN_REF} -dbtype nucl -parse_seqids

#ml LAST
#lastdb -uNEAR -R01 ${IN_REF} ${IN_REF}
