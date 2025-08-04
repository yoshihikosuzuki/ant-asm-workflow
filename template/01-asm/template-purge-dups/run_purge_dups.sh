#!/bin/bash
#SBATCH -J purge_dups
#SBATCH -o purge_dups.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=100G
#SBATCH -t 24:00:00
source ../../config.sh
set -eu
module load ${_PURGE_DUPS} ${_SEQKIT}
set -x

# According to `calcuts` command:
#     -l    INT      lower bound for read depth
#     -m    INT      transition between haploid and diploid
#     -u    INT      upper bound for read depth
L=
M=
U=

if [ -z "${L}" ]; then
    L=$(cut -f1,1 cutoffs)
fi
if [ -z "${M}" ]; then
    M=$(cut -f4,4 cutoffs)
fi
if [ -z "${U}" ]; then
    U=$(cut -f6,6 cutoffs)
fi

CONTIGS=contigs.fasta
PAF_CONTIGS=contigs.paf
CUTOFFS=cutoffs

calcuts -l ${L} -m ${M} -u ${U} PB.stat >${CUTOFFS}
purge_dups -2 -T ${CUTOFFS} -c PB.base.cov ${PAF_CONTIGS} >dups.bed
get_seqs dups.bed ${CONTIGS}
echo "Finished purge_dups"
seqkit stats -a purged.fa >purged.fa.stats
echo "Purged contig stats:"
cat purged.fa.stats
