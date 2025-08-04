#!/bin/bash
#SBATCH -J curation
#SBATCH -o curation.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=500G
#SBATCH -t 10:00:00
source ../../config.sh
set -eu
module load ${_SAMTOOLS} ${_BWA} ${_3DDNA}
set -x

SCAFS=scaffolds.fasta
REVIEW_ASSEMBLY=scaffolds.review.assembly
MERGED_NODUPS=merged_nodups.txt

3d-dna-post-review -r ${REVIEW_ASSEMBLY} ${SCAFS} ${MERGED_NODUPS}
