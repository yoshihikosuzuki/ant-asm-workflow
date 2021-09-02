#!/bin/bash
#SBATCH -J peregrine
#SBATCH -o peregrine.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 128
#SBATCH --mem=500G
#SBATCH -t 48:00:00
shopt -s expand_aliases && source ~/.bashrc && set -e || exit 1
source ../../config.sh

IN_FASTX=hifi.fastq
N_THREAD=128
OUT_DIR=asm

N_THREAD_SETTING="${N_THREAD} ${N_THREAD} ${N_THREAD} ${N_THREAD} ${N_THREAD} ${N_THREAD} ${N_THREAD} ${N_THREAD} ${N_THREAD}"

ml Other/peregrine Other/seqkit

find $PWD -name "${IN_FASTX}" > reads.fofn
peregrine reads.fofn ${N_THREAD_SETTING} --with-consensus --shimmer-r 3 --best_n_ovlp 8 --output ${OUT_DIR}
echo "Contig stats (${OUT_DIR}/p_ctg_cns.fa):"
seqkit stats -a ${OUT_DIR}/p_ctg_cns.fa

if [ "$AUTO_DEL" = "true" ]; then
    source ./remove_tmp_files.sh
fi
