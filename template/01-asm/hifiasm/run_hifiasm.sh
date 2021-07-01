#!/bin/bash
#SBATCH -J hifiasm
#SBATCH -o hifiasm.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 128
#SBATCH --mem=500G
#SBATCH -t 24:00:00

IN_FASTX=
HIC_READS_1=
HIC_READS_2=
N_THREADS=128

OUT_PREFIX=${IN_FASTX%.gz}
OUT_PREFIX=${OUT_PREFIX%.*}

ml hifiasm
## Normal mode
hifiasm -o ${OUT_PREFIX} -t ${N_THREADS} ${IN_FASTX}
## Hi-C mode
hifiasm -o ${OUT_PREFIX} -t ${N_THREADS} --h1 ${HIC_READS_1} --h2 ${HIC_READS_2} ${IN_FASTX}

ml gfatools
for DATA in *tg.gfa; do
    gfatools gfa2fa ${DATA} > ${DATA%.gfa}.fasta
done
for DATA in *.noseq.gfa; do
    awk 'BEGIN {print "Name,Depth"} /^S/ {print $2 "," substr($5,6)}' ${DATA} > ${DATA%.noseq.gfa}.depth.csv
done
