#!/bin/bash
#SBATCH -J deepvariant
#SBATCH -o deepvariant.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 128
#SBATCH --mem=500G
#SBATCH -t 24:00:00
shopt -s expand_aliases && source ~/.bashrc && set -e || exit 1

REF=scaffolds.fasta
READS=hifi.fastq
N_THREADS=128

# NOTE: Assuming the specific directory structure for input BAM file
_REF=$(basename ${REF} .gz)
_READS=$(basename ${READS} .gz)
OUT_PREFIX=${_REF%.*}.${_READS%.*}.winnowmap
OUT_BAM=${OUT_PREFIX}.sorted.bam
IN_SORTED_BAM=../04-winnowmap/${OUT_BAM}

OUT_PREFIX=${_REF%.*}.${_READS%.*}.deepvariant
OUT_VCF=${OUT_PREFIX}.vcf

#ml samtools
#samtools faidx ${REF}
#samtools index -@${N_THREADS} ${IN_SORTED_BAM}

ml deepvariant

run_deepvariant --num_shards ${N_THREADS} --model_type PACBIO --ref ${REF} --reads ${IN_SORTED_BAM} --output_vcf ${OUT_VCF}
