#!/bin/bash
#SBATCH -J peregrine
#SBATCH -o peregrine.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 128
#SBATCH --mem=500G
#SBATCH -t 10:00:00
shopt -s expand_aliases && source ~/.bashrc && set -e || exit 1

IN_FASTX=
N_THREAD=128

N_THREAD_SETTING="${N_THREAD} ${N_THREAD} ${N_THREAD} ${N_THREAD} ${N_THREAD} ${N_THREAD} ${N_THREAD} ${N_THREAD} ${N_THREAD}"

ml peregrine
find $PWD -name "${IN_FASTX}" > reads.fofn
run_peregrine reads.fofn ${N_THREAD_SETTING} --with-consensus --shimmer-r 3 --best_n_ovlp 8 --output asm
