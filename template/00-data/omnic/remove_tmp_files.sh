#!/bin/bash
#SBATCH -J remove_tmp_files
#SBATCH -o remove_tmp_files.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=1G
#SBATCH -t 1:00:00
shopt -s expand_aliases && source ~/.bashrc && set -e || exit 1
source ../../config.sh

FASTK_PREFIX=omnic.fastk

ml ${_FASTK}

Fastrm ${FASTK_PREFIX}
rm -rf tmp/
