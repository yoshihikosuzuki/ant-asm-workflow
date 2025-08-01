#!/bin/bash
#SBATCH -J omnic-qc
#SBATCH -o run_all.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=1G
#SBATCH -t 72:00:00
source ../../config.sh
set -eu
ml ${_FASTK}
set -x

FASTK=$(sbatch run_fastk.sh | cut -f 4 -d' ')
GENESCOPE=$(sbatch -d afterany:${FASTK} run_genescope.sh | cut -f 4 -d' ')
WAIT_JOBS=${GENESCOPE}

if [ "$AUTO_DEL" = "true" ]; then
    WAIT_JOBS=$(sbatch -d afterany:${WAIT_JOBS} remove_tmp_files.sh | cut -f 4 -d' ')
fi

srun -p compute -c 1 --mem 1G -t 1:00:00 -d afterany:${WAIT_JOBS} --wait=0 sleep 1s
