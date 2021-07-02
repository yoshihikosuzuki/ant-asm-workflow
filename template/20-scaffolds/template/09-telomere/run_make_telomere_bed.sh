#!/bin/bash
#SBATCH -J telomere
#SBATCH -o telomere.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=100G
#SBATCH -t 24:00:00

IN_FASTA=scaffolds.fasta
TELOMERE_MOTIF="TTAGGG"

ml TRF #make_telomere_bed

make_telomere_bed ${IN_FASTA} ${TELOMERE_MOTIF}
