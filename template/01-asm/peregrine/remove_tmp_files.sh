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

find asm/ -maxdepth 1 -name "[0-3]*" | while read DIR; do rm -rf ${DIR}; done
find asm/4-cns/ -name "cns-[0-9]*" | while read DIR; do rm -rf ${DIR}; done
rm -rf asm/4-cns/ctg_index/ asm/4-cns/map-* asm/4-cns/tmp/
