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

PAF_READS=reads.paf

ml purge_dups

pbcstat ${PAF_READS}
hist_plot.py PB.stat PB.hist.png
calcuts PB.stat > cutoffs
info echo -n "Finished calcuts: "
cat cutoffs
