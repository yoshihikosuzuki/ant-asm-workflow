#!/bin/bash
#SBATCH -J hicanu
#SBATCH -o hicanu.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 64
#SBATCH --mem=500G
#SBATCH -t 24:00:00

IN_FASTX=
OUT_PREFIX=
GENOME_SIZE=

ml canu

canu -d asm -p ${OUT_PREFIX} genomeSize=${GENOME_SIZE} useGrid=false -pacbio-hifi ${IN_FASTX}
