#!/bin/bash
#SBATCH -J hicanu
#SBATCH -o hicanu.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 128
#SBATCH --mem=500G
#SBATCH -t 48:00:00
shopt -s expand_aliases && source ~/.bashrc && set -e || exit 1
source ../../config.sh

IN_FASTX=hifi.fastq

OUT_PREFIX=$(basename ${IN_FASTX} .gz)
OUT_PREFIX=${OUT_PREFIX%.*}.hicanu

ml ${_CANU} ${_SEQKIT}

canu -d . -p ${OUT_PREFIX} genomeSize=${GENOME_SIZE} useGrid=false -pacbio-hifi ${IN_FASTX}
echo "Contig stats (${OUT_PREFIX}.contigs.fasta):"
seqkit stats -a ${OUT_PREFIX}.contigs.fasta

if [ "$AUTO_DEL" = "true" ]; then
    source ./remove_tmp_files.sh
fi
