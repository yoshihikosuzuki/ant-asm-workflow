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

rm -rf tmp/
find RUN/ -maxdepth 1 -name "0[0-9]*" | while read DIR; do rm -rf $DIR; done
find RUN/ -maxdepth 1 -name "1[0-8]*" | while read DIR; do rm -rf $DIR; done

