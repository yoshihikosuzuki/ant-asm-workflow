#!/bin/bash
#SBATCH -J hifi-qc
#SBATCH -o run_all.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=40G
#SBATCH -t 72:00:00
shopt -s expand_aliases && source ~/.bashrc && set -e || exit 1

FASTK=$(sbatch run_fastk.sh | cut -f 4 -d' ')
GENESCOPE=$(sbatch -d afterany:${FASTK} run_genescope.sh | cut -f 4 -d' ')
GENOMESCOPE=$(sbatch run_genomescope.sh | cut -f 4 -d' ')
srun -p compute -c 1 --mem 1G -t 1:00:00 -d afterany:${GENESCOPE},${GENOMESCOPE} --wait=0 sleep 1s