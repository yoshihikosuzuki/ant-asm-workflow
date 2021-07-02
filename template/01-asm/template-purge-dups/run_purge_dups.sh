#!/bin/bash
#SBATCH -J purge_dups
#SBATCH -o purge_dups.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=100G
#SBATCH -t 24:00:00
shopt -s expand_aliases && source ~/.bashrc && set -e || exit 1

L=
M=
U=

PAF_CONTIGS=contigs.paf
CUTOFFS=cutoffs

ml purge_dups

calcuts -l ${L} -m ${M} -u ${U} PB.stat > ${CUTOFFS}
purge_dups -2 -T ${CUTOFFS} -c PB.base.cov ${PAF_CONTIGS} > dups.bed
get_seqs dups.bed ${CONTIGS}
info echo "Finished purge_dups"
