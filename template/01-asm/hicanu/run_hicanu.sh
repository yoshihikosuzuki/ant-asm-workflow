#!/bin/bash
#SBATCH -J hicanu
#SBATCH -o hicanu.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 128
#SBATCH --mem=500G
#SBATCH -t 48:00:00
shopt -s expand_aliases && source ~/.bashrc && set -e || exit 1

IN_FASTX=hifi.fastq
GENOME_SIZE=300000000

OUT_PREFIX=${IN_FASTX%.gz}
OUT_PREFIX=${OUT_PREFIX%.*}.hicanu

ml Other/canu

canu -d . -p ${OUT_PREFIX} genomeSize=${GENOME_SIZE} useGrid=false -pacbio-hifi ${IN_FASTX}
