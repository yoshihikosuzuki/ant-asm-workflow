#!/bin/bash
#SBATCH -J mummer
#SBATCH -o mummer.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 128
#SBATCH --mem=100G
#SBATCH -t 10:00:00
shopt -s expand_aliases && source ~/.bashrc && set -e || exit 1
source ../../../config.sh

REF=contigs.hap1.fasta
QUERY=contigs.hap2.fasta
N_THREADS=128

_REF=$(basename ${REF} .gz)
_QUERY=$(basename ${QUERY} .gz)
OUT_PREFIX=${_REF%.*}.${_QUERY%.*}.mummer
OUT_DELTA=${OUT_PREFIX}.delta
OUT_FILT_DELTA=${OUT_PREFIX}.filtered.delta
OUT_GP=${OUT_PREFIX}.gp

ml ${_MUMMER}

set term eps
nucmer -t ${N_THREADS} -p ${OUT_PREFIX} ${REF} ${QUERY}
## Only 1-to-1 alignments
delta-filter -1 -i 90 -l 1000 ${OUT_DELTA} >${OUT_FILT_DELTA}
## Keep secondary alignments
#delta-filter -i 90 -l 1000 ${OUT_DELTA} >${OUT_FILT_DELTA}
mummerplot -R ${REF} -Q ${QUERY} -p ${OUT_PREFIX} --postscript --layout --fat ${OUT_FILT_DELTA}

# Change font and font size
sed -i 's/"Courier" 8/"Helvetica" 5/' ${OUT_GP}
# Write characters of labels as they are
# TODO: Change `1,1` to e.g. `0.7,1` if aspect ratio is wrong
sed -i 's/set size 1,1/set termoption noenhanced\nset size 1,1/' ${OUT_GP}
# Use solid lines for contig boundaries
sed -i 's/set grid/set grid linetype 1 linewidth 0.1 linecolor "light-gray"/' ${OUT_GP}
# Remove contig boundaries
#sed -i 's/set grid/#set grid/' ${OUT_GP}
# Remove contig names
#sed -i 's/unset key/unset key\nunset xtics\nunset ytics/' ${OUT_GP}
# Change axis titles
sed -i "s/set xlabel \"REF\"/set xlabel \"${REF}\"/" ${OUT_GP}
sed -i "s/set ylabel \"QRY\"/set ylabel \"${QUERY}\"/" ${OUT_GP}
# Use thinner dots and lines for alignments (lw: linewidth，pt: pointtype，ps: pointsize)
sed -i 's/lw 2 pt 6 ps 0.5/lw 0.5 pt 7 ps 0.05/g' ${OUT_GP}
gnuplot ${OUT_GP}
