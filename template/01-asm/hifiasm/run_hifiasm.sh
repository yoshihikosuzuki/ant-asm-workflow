#!/bin/bash
#SBATCH -J hifiasm
#SBATCH -o hifiasm.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 128
#SBATCH --mem=500G
#SBATCH -t 24:00:00
source ../../config.sh
set -eu
module load ${_HIFIASM} ${_GFATOOLS} ${_SEQKIT}
set -x

IN_FASTX=hifi.fastq${HIFI_GZ}
N_THREADS=128

OUT_PREFIX=$(basename ${IN_FASTX} .gz)
OUT_PREFIX=${OUT_PREFIX%.*}.hifiasm

hifiasm -o ${OUT_PREFIX} -t ${N_THREADS} ${IN_FASTX}

for ASM in *tg.gfa; do
    gfatools gfa2fa ${ASM} >${ASM%.gfa}.fasta
done

for ASM in *.bp.p_ctg.fasta; do
    seqkit stats -a ${ASM} >${ASM}.stats
    echo "Contig stats (${ASM}.stats):"
    cat ${ASM}.stats
done

if [ "$AUTO_DEL" = "true" ]; then
    source ./remove_tmp_files.sh
fi
