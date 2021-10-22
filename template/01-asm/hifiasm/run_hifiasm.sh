#!/bin/bash
#SBATCH -J hifiasm
#SBATCH -o hifiasm.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 128
#SBATCH --mem=500G
#SBATCH -t 24:00:00
shopt -s expand_aliases && source ~/.bashrc && set -e || exit 1
source ../../config.sh

IN_FASTX=hifi.fastq
N_THREADS=128

OUT_PREFIX=$(basename ${IN_FASTX} .gz)
OUT_PREFIX=${OUT_PREFIX%.*}.hifiasm

ml ${_HIFIASM} ${_GFATOOLS} ${_SEQKIT}

hifiasm -o ${OUT_PREFIX} -t ${N_THREADS} ${IN_FASTX}
for DATA in *tg.gfa; do
    gfatools gfa2fa ${DATA} > ${DATA%.gfa}.fasta
done

P_UTG=${OUT_PREFIX}.bp.p_utg.fasta
P_CTG=${OUT_PREFIX}.bp.p_ctg.fasta
for DATA in ${P_UTG} ${P_CTG}; do
    echo "Contig stats (${DATA}):"
    seqkit stats -a ${DATA}
done

if [ "$AUTO_DEL" = "true" ]; then
    source ./remove_tmp_files.sh
fi
