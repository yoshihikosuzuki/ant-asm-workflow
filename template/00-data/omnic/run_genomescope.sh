#!/bin/bash
#SBATCH -J genomescope
#SBATCH -o genomescope.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=100G
#SBATCH -t 48:00:00
shopt -s expand_aliases && source ~/.bashrc && set -e || exit 1
source ../../config.sh

IN_FNAMES="omnic_R1_001.fastq omnic_R2_001.fastq"
IN_FOFN=omnic.fofn
find ${IN_FNAMES} > ${IN_FOFN}

K=21
PLOIDY=2
HIST_MAX=10000
N_THREADS=16
GB_MEMORY=30
TMP_DIR=tmp

OUT_PREFIX=${IN_FOFN/.fofn/.genomescope}
OUT_HIST=${OUT_PREFIX}.hist

ml Other/genomescope Other/smudgeplot

mkdir -p ${TMP_DIR}
kmc -k${K} -m${GB_MEMORY} -ci1 -cs${HIST_MAX} -t${N_THREADS} -fm @${IN_FOFN} ${OUT_PREFIX} ${TMP_DIR}
kmc_tools transform ${OUT_PREFIX} -cx${HIST_MAX} histogram ${OUT_HIST}
genomescope.R -i ${OUT_HIST} -o ${OUT_PREFIX} -p${PLOIDY} -k ${K}

L=$(smudgeplot.py cutoff ${OUT_HIST} L)
U=$(smudgeplot.py cutoff ${OUT_HIST} U)
echo L=$L U=$U
# NOTE: L is typically within 20-200 and U is 500-3000
OUT_KMERS=${OUT_PREFIX}_L${L}_U${U}
kmc_tools transform ${OUT_PREFIX} -ci${L} -cx${U} reduce ${OUT_KMERS}
smudge_pairs ${OUT_KMERS} ${OUT_KMERS}_coverages.tsv ${OUT_KMERS}_pairs.tsv >${OUT_KMERS}_familysizes.tsv
smudgeplot.py plot -o ${OUT_KMERS} ${OUT_KMERS}_coverages.tsv
# NOTE: To explicitly specify the haploid depth, run below instead
#HAPLO_DEPTH=
#smudgeplot.py plot -o ${OUT_KMERS} -n ${HAPLO_DEPTH} ${OUT_KMERS}_coverages.tsv
