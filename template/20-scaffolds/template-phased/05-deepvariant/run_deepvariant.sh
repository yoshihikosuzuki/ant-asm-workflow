#!/bin/bash
#SBATCH -J deepvariant
#SBATCH -o deepvariant.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 128
#SBATCH --mem=500G
#SBATCH -t 48:00:00
shopt -s expand_aliases && source ~/.bashrc && set -e || exit 1
source ../../../config.sh

REF=scaffolds.fasta
READS=hifi.fastq
IN_SORTED_BAM=scaffolds.hifi.winnowmap.sorted.bam
N_THREADS=128

OUT_PREFIX=${REF%.*}.${READS%.*}.deepvariant
OUT_VCF=${OUT_PREFIX}.vcf

ml ${_DEEPVARIANT}

run_deepvariant --num_shards ${N_THREADS} --model_type PACBIO --ref ${REF} --reads ${IN_SORTED_BAM} --output_vcf ${OUT_VCF}
