#!/bin/bash
#SBATCH -J omnic-qc
#SBATCH -o run_all.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=40G
#SBATCH -t 72:00:00
shopt -s expand_aliases && source ~/.bashrc && set -e || exit 1
source ../../config.sh

FASTK=$(sbatch run_fastk.sh | cut -f 4 -d' ')
GENESCOPE=$(sbatch -d afterany:${FASTK} run_genescope.sh | cut -f 4 -d' ')
GENOMESCOPE=$(sbatch run_genomescope.sh | cut -f 4 -d' ')
if [ "$AUTO_DEL" = "true" ]; then
    DEL=,$(sbatch -d afterany:${GENESCOPE},${GENOMESCOPE} remove_tmp_files.sh | cut -f 4 -d' ')
fi
srun -p compute -c 1 --mem 1G -t 1:00:00 -d afterany:${GENESCOPE},${GENOMESCOPE}${DEL} --wait=0 sleep 1s
