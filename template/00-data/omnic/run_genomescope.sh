#!/bin/bash
#SBATCH -J genomescope
#SBATCH -o genomescope.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=100G
#SBATCH -t 10:00:00

IN_FASTA=

_OUT_PREFIX=${IN_FASTA%.gz}
OUT_PREFIX=${_OUT_PREFIX%.*}.genomescope
OUT_HIST=${OUT_PREFIX}.hist

# NOTE: K=21 and HIST_MAX=10000 are recommended values by GenomeScope
K=21
HIST_MAX=10000
N_THREADS=16
GB_MEMORY=30

ml genomescope

mkdir -p tmp
kmc -k${K} -m${GB_MEMORY} -ci1 -cs${HIST_MAX} -t${N_THREADS} -fm ${IN_FASTA} ${OUT_PREFIX} tmp
kmc_tools transform ${OUT_PREFIX} -cx${HIST_MAX} histogram ${OUT_HIST}
genomescope.R -i ${OUT_HIST} -o ${OUT_PREFIX} -k ${K}

ml smudgeplot

L=$(smudgeplot.py cutoff ${OUT_HIST} L)
U=$(smudgeplot.py cutoff ${OUT_HIST} U)
echo L=$L U=$U
# L should be like 20 - 200
# U should be like 500 - 3000
OUT_KMERS=${OUT_PREFIX}_L${L}_U${U}
kmc_tools transform ${OUT_PREFIX} -ci${L} -cx${U} reduce ${OUT_KMERS}
smudge_pairs ${OUT_KMERS} ${OUT_KMERS}_coverages.tsv ${OUT_KMERS}_pairs.tsv >${OUT_KMERS}_familysizes.tsv
smudgeplot.py plot -o ${OUT_KMERS} ${OUT_KMERS}_coverages.tsv
#HAPLO_DEPTH=
#smudgeplot.py plot -o ${OUT_KMERS} -n ${HAPLO_DEPTH} ${OUT_KMERS}_coverages.tsv
