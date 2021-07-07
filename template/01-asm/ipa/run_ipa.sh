#!/bin/bash
#SBATCH -J ipa
#SBATCH -o ipa.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 128
#SBATCH --mem=500G
#SBATCH -t 48:00:00
shopt -s expand_aliases && source ~/.bashrc && set -e || exit 1

IN_FASTX=hifi.fastq
N_THREAD=128
TMP_DIR=tmp

ml Other/pbipa

mkdir -p ${TMP_DIR}
ipa local --input-fn ${IN_FASTX} --nthreads ${N_THREAD} --tmp-dir ${TMP_DIR} --verbose
