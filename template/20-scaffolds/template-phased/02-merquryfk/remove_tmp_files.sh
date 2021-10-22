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
source ../../../config.sh

ml ${_FASTK}

Fastrm hifi.fastk scaffolds.hap1.fastk scaffolds.hap2.fastk scaffolds.hap1.fastk.hifi.fastk scaffolds.hap2.fastk.hifi.fastk
rm -rf tmp/