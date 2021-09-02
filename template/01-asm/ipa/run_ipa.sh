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
source ../../config.sh

IN_FASTX=hifi.fastq
N_THREAD=128
TMP_DIR=$(readlink -f tmp/)

ml Other/pbipa Other/seqkit

mkdir -p ${TMP_DIR}
ipa local --input-fn ${IN_FASTX} --nthreads ${N_THREAD} --tmp-dir ${TMP_DIR} --verbose
echo "Contig stats (RUN/assembly-results/final.p_ctg.fasta):"
seqkit stats -a RUN/assembly-results/final.p_ctg.fasta

if [ "$AUTO_DEL" = "true" ]; then
    source ./remove_tmp_files.sh
fi
