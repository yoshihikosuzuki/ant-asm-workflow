#!/bin/bash
#SBATCH -J make_index
#SBATCH -o make_index.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=100G
#SBATCH -t 10:00:00
shopt -s expand_aliases && source ~/.bashrc && set -e || exit 1
source ../../config.sh

IN_REF=contigs.fasta

ml samtools bwa

samtools faidx ${IN_REF}
bwa index -p ${IN_REF} ${IN_REF}
