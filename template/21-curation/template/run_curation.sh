#!/bin/bash
#SBATCH -J curation
#SBATCH -o curation.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=500G
#SBATCH -t 10:00:00
shopt -s expand_aliases && source ~/.bashrc && set -e || exit 1
source ../../config.sh

SCAFS=draft.fasta
REVIEW_ASSEMBLY=scaffolds.review.assembly
MERGED_NODUPS=merged_nodups.txt

ml ${_SAMTOOLS} ${_BWA} ${_3DDNA}

3d-dna-post-review -r ${REVIEW_ASSEMBLY} ${SCAFS} ${MERGED_NODUPS}
