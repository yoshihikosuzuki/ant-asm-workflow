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

ml Other/hifiasm Other/gfatools Other/seqkit

hifiasm -o ${OUT_PREFIX} -t ${N_THREADS} --h1 ${HIC_READS_1} --h2 ${HIC_READS_2} ${IN_FASTX}
for DATA in *tg.gfa; do
    gfatools gfa2fa ${DATA} > ${DATA%.gfa}.fasta
done

echo "Contig stats (${OUT_PREFIX}.bp.p_utg.fasta):"
seqkit stats -a ${OUT_PREFIX}.bp.p_utg.fasta
echo "Contig stats (${OUT_PREFIX}.bp.p_ctg.fasta):"
seqkit stats -a ${OUT_PREFIX}.bp.p_ctg.fasta

if [ "$AUTO_DEL" = "true" ]; then
    source ./remove_tmp_files.sh
fi
