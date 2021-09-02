#!/bin/bash
#SBATCH -J purge_dups_plot
#SBATCH -o purge_dups_plot.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 128
#SBATCH --mem=500G
#SBATCH -t 24:00:00
shopt -s expand_aliases && source ~/.bashrc && set -e || exit 1
source ../../config.sh

CONTIGS=contigs.fasta
READS=hifi.fastq
N_THREADS=128

CONTIGS_PREFIX=${CONTIGS%.gz}
CONTIGS_SPLIT=${CONTIGS_PREFIX%.*}.split.fa

ml Other/minimap2 Other/purge_dups

minimap2 -t${N_THREADS} -xmap-pb ${CONTIGS} ${READS} > reads.paf
split_fa ${CONTIGS} > ${CONTIGS_SPLIT}
minimap2 -t${N_THREADS} -xasm5 -DP ${CONTIGS_SPLIT} ${CONTIGS_SPLIT} > contigs.paf
pbcstat reads.paf
hist_plot.py PB.stat PB.hist.png
calcuts PB.stat > cutoffs
echo -n "Automatically estimated cutoff values: "
cat cutoffs

if [ "$AUTO_DEL" = "true" ]; then
    source ./remove_tmp_files.sh
fi
