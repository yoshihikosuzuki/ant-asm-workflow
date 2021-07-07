#!/bin/bash
#SBATCH -J telomere
#SBATCH -o telomere.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=100G
#SBATCH -t 24:00:00

IN_FASTA=contigs.fasta
TELOMERE_MOTIF="TTAGG"

ml Other/make_telomere_bed

make_telomere_bed ${IN_FASTA} ${TELOMERE_MOTIF}

# NOTE: If TRF freezes, run the followings instead:
#rm -f *.trf
#make_telomere_bed -s ${IN_FASTA} ${TELOMERE_MOTIF}
