#!/bin/bash
#SBATCH -J reference
#SBATCH -o reference.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=500G
#SBATCH -t 24:00:00
shopt -s expand_aliases && source ~/.bashrc && set -e || exit 1

IN_REF=

ml samtools
samtools faidx ${IN_REF}

ml bwa
bwa index -p ${IN_REF}.bwa ${IN_REF}

ml blast+
makeblastdb -in ${IN_REF} -out ${IN_REF}.blast -dbtype nucl -parse_seqids

ml LAST
lastdb -uNEAR -cR01 ${IN_REF}.last ${IN_REF}

ml pbmm2
IN_FASTA=${IN_REF%.*}.fasta
if [ ${IN_REF##*.} != "fasta" ]; then
    ln -sf ${IN_REF} ${IN_FASTA}
fi
pbmm2 index ${IN_FASTA} ${IN_FASTA}.mmi
