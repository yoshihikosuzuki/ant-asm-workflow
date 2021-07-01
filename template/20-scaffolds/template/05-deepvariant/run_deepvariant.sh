#!/bin/bash
#SBATCH -J deepvariant
#SBATCH -o deepvariant.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 128
#SBATCH --mem=500G
#SBATCH -t 10:00:00
shopt -s expand_aliases && source ~/.bashrc && set -e || exit 1

REF=
IN_SORTED_BAM=
OUT_VCF=
N_THREADS=128

ml samtools
samtools faidx ${REF}
samtools index -@${N_THREADS} ${IN_SORTED_BAM}

ml deepvariant
run_deepvariant --num_shards ${N_THREADS} --model_type PACBIO --ref ${REF} --reads ${IN_SORTED_BAM} --output_vcf ${OUT_VCF}
