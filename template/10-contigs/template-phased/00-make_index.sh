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

ml ${_SAMTOOLS} ${_BWA} ${_3DDNA}

for IN_REF in contigs.fasta contigs.hap1.fasta contigs.hap2.fasta; do
    _PREFIX=${IN_REF%.gz}
    _PREFIX=${_PREFIX%.*}
    OUT_ASSEMBLY=${_PREFIX}.assembly
    OUT_CHROM_SIZES=${_PREFIX}.chrom_sizes
    samtools faidx ${IN_REF}
    bwa index -p ${IN_REF} ${IN_REF}
    3d-dna-fasta2assembly ${IN_REF} >${OUT_ASSEMBLY}
    awk 'NF == 3 {print substr($1,2) "\t" $3}' ${OUT_ASSEMBLY} >${OUT_CHROM_SIZES}
done
