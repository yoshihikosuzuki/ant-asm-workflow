#!/bin/bash
#SBATCH -J scaf-eval
#SBATCH -o run_all.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=1G
#SBATCH -t 72:00:00
shopt -s expand_aliases && source ~/.bashrc && set -e || exit 1

HAP1=scaffolds.hap1.fasta
HAP2=scaffolds.hap2.fasta
MERGED=scaffolds.fasta

cat ${HAP1} ${HAP2} >${MERGED}

MAKE_INDEX=$(sbatch 00-make_index.sh | cut -f 4 -d' ')

cd 01-busco-all/ &&
    BUSCO_ALL=$(sbatch run_busco.sh | cut -f 4 -d' ') &&
    cd ..
cd 01-busco-hap1/ &&
    BUSCO1=$(sbatch run_busco.sh | cut -f 4 -d' ') &&
    cd ..
cd 01-busco-hap2/ &&
    BUSCO2=$(sbatch run_busco.sh | cut -f 4 -d' ') &&
    cd ..
cd 02-merquryfk/ &&
    MERQURYFK=$(sbatch run_merquryfk.sh | cut -f 4 -d' ') &&
    cd ..
cd 04-winnowmap/ &&
    WINNOWMAP=$(sbatch -d afterany:${MAKE_INDEX} run_winnowmap.sh | cut -f 4 -d' ') &&
    cd ..
cd 05-deepvariant/ &&
    DEEPVARIANT=$(sbatch -d afterany:${WINNOWMAP} run_deepvariant.sh | cut -f 4 -d' ') &&
    cd ..
cd 06-mapqv/ &&
    MAPQV=$(sbatch -d afterany:${DEEPVARIANT} run_mapqv.sh | cut -f 4 -d' ') &&
    cd ..
cd 07-asset/ &&
    ASSET=$(sbatch -d afterany:${WINNOWMAP} run_asset.sh | cut -f 4 -d' ') &&
    cd ..
cd 09-telomere/ &&
    TELOMERE=$(sbatch run_make_telomere_bed.sh | cut -f 4 -d' ') &&
    cd ..
WAIT_JOBS=${BUSCO_ALL},${BUSCO1},${BUSCO2},${MERQURYFK},${MAPQV},${ASSET},${TELOMERE}

srun -p compute -c 1 --mem 1G -t 1:00:00 -d afterany:${WAIT_JOBS} --wait=0 sleep 1s
