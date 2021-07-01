#!/bin/bash
#SBATCH -J telomere
#SBATCH -o telomere.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=100G
#SBATCH -t 24:00:00

TELOMERE_MOTIF="TTAGGG"

IN_FASTAS=(   # <fasta_file> per line

)

ml TRF

for IN_FASTA in "${IN_FASTAS[@]}"; do
    make_telomere_bed ${IN_FASTA} ${TELOMERE_MOTIF}
done
