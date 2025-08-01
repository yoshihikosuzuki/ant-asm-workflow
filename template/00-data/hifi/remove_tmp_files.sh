#!/bin/bash
#SBATCH -J remove_tmp_files
#SBATCH -o remove_tmp_files.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=1G
#SBATCH -t 1:00:00
source ../../config.sh
set -eu
ml ${_FASTK}
set -x

FASTK_PREFIX=hifi.fastk

Fastrm ${FASTK_PREFIX}
rm -rf tmp/
