#!/bin/bash
#SBATCH -J ipa
#SBATCH -o ipa.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 128
#SBATCH --mem=500G
#SBATCH -t 24:00:00
shopt -s expand_aliases && source ~/.bashrc && set -e || exit 1

IN_FASTX=
N_THREAD=128

ml pbipa
mkdir -p tmp
ipa local --input-fn ${IN_FASTX} --nthreads ${N_THREAD} --tmp-dir tmp --verbose
