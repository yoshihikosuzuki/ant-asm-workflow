#!/bin/bash
#SBATCH -J telomere
#SBATCH -o telomere.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=100G
#SBATCH -t 24:00:00
shopt -s expand_aliases && source ~/.bashrc && set -e || exit 1
source ../../../config.sh

IN_FASTA=scaffolds.fasta
MIN_NCOPY=${TELOMERE_MIN_NCOPY}

OUT_BED=${IN_FASTA/.fasta/.telomere.bed}
OUT_FILT_BED=${OUT_BED/.bed/.filtered.bed}

ml ${_MAKE_TELOMERE_BED}

make_telomere_bed ${IN_FASTA} ${TELOMERE_MOTIF}
# NOTE: If TRF freezes, run the followings instead:
#rm -f *.trf
#make_telomere_bed -s ${IN_FASTA} ${TELOMERE_MOTIF}


filter_bed -m ${MIN_NCOPY} ${OUT_BED} > ${OUT_FILT_BED}

# Generate .bed files for JBAT

bed_to_jbat() {
    CHROM_SIZES=$1
    BED_FILE=$2
    awk 'BEGIN {l=0} FNR == NR {offset[$1]=l; l+=$2; next} {printf "assembly\t" offset[$1]+$2 "\t" offset[$1]+$3; for(i=4;i<=NF;i++) printf "\t" $i; print ""}' ${CHROM_SIZES} ${BED_FILE}
}

bed_to_jbat ../${IN_FASTA/.fasta/.chrom_sizes} ${OUT_BED} > ${OUT_BED/.bed/.JBAT.bed}
bed_to_jbat ../${IN_FASTA/.fasta/.chrom_sizes} ${OUT_FILT_BED} > ${OUT_FILT_BED/.bed/.JBAT.bed}
