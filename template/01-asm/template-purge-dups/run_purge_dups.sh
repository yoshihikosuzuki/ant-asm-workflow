#!/bin/bash
#SBATCH -J purge_dups
#SBATCH -o purge_dups.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 128
#SBATCH --mem=500G
#SBATCH -t 24:00:00
shopt -s expand_aliases && source ~/.bashrc && set -e || exit 1

CONTIGS=
READS=
OUT_PREFIX=
# EXAMPLE: 6G
GENOME_SIZE=
N_THREADS=128

ml minimap2 purge_dups

calc_cutoffs () {
    minimap2 -t${N_THREADS} -I${GENOME_SIZE} -xmap-pb ${CONTIGS} ${READS} > ${OUT_PREFIX}.reads.paf
    split_fa ${CONTIGS} > ${OUT_PREFIX}.split.fa
    minimap2 -t${N_THREADS} -I${GENOME_SIZE} -xasm5 -DP ${OUT_PREFIX}.split.fa ${OUT_PREFIX}.split.fa > ${OUT_PREFIX}.contigs.paf
    pbcstat ${OUT_PREFIX}.reads.paf
    hist_plot.py PB.stat PB.hist.png
    calcuts PB.stat > ${OUT_PREFIX}.cutoffs
    info echo -n "Finished calcuts: "
    cat ${OUT_PREFIX}.cutoffs
}

run_purge () {   # Usage: `$ run_purge <EH threshold> <HD threshold> <DR threshold>`
    EH=$1
    HD=$2
    DR=$3
    calcuts -l $EH -m $HD -u $DR PB.stat > ${OUT_PREFIX}.cutoffs
    purge_dups -2 -T ${OUT_PREFIX}.cutoffs -c PB.base.cov ${OUT_PREFIX}.contigs.paf > ${OUT_PREFIX}.dups.bed
    get_seqs ${OUT_PREFIX}.dups.bed ${CONTIGS}
    info echo "Finished purge_dups"
}

calc_cutoffs
# NOTE: see the `PB.hist.png` and change the thresholds
#run_purge 5 20 100
