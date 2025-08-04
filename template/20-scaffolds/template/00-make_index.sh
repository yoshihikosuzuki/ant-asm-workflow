#!/bin/bash
#SBATCH -J make_index
#SBATCH -o make_index.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=100G
#SBATCH -t 10:00:00
source ../../config.sh
set -eu
module load ${_SAMTOOLS} ${_BWA} ${_3DDNA}
set -x

IN_REF=scaffolds.fasta
OUT_ASSEMBLY=${IN_REF%.*}.assembly
OUT_CHROM_SIZES=${IN_REF%.*}.chrom_sizes

samtools faidx ${IN_REF}
bwa index -p ${IN_REF} ${IN_REF}
3d-dna-fasta2assembly ${IN_REF} >${OUT_ASSEMBLY}
awk 'NF == 3 {print substr($1,2) "\t" $3}' ${OUT_ASSEMBLY} >${OUT_CHROM_SIZES}
