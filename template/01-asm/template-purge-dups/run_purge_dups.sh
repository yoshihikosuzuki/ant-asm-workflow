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
source ../../config.sh

# According to `calcuts` command:
#     -l    INT      lower bound for read depth
#     -m    INT      transition between haploid and diploid
#     -u    INT      upper bound for read depth
L=
M=
U=

CONTIGS=contigs.fasta
PAF_CONTIGS=contigs.paf
CUTOFFS=cutoffs

ml ${_PURGE_DUPS} ${_SEQKIT}

calcuts -l ${L} -m ${M} -u ${U} PB.stat > ${CUTOFFS}
purge_dups -2 -T ${CUTOFFS} -c PB.base.cov ${PAF_CONTIGS} > dups.bed
get_seqs dups.bed ${CONTIGS}
echo "Finished purge_dups"
echo "Purged contig stats:"
seqkit stats -a purged.fa
