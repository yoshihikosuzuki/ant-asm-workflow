#!/bin/bash
set -eux

PATH_TO_TEMPLATE=./template
GZ=".gz"
# GZ=""

###################################################################################################
## make directory and download test data
###################################################################################################
rm -rf test
cp -r ${PATH_TO_TEMPLATE} test
cd test
wget -O - https://mlab.cb.k.u-tokyo.ac.jp/~yoshihiko_s/ant-asm-workflow/reads.tar.gz | tar xzvf -

###################################################################################################
## `00-data`: make symlink to reads of the test data
###################################################################################################
cd 00-data

cd hifi
ln -sf ../../reads/hifi.fastq${GZ} .
sbatch run_all.sh
cd ..

cd omnic
ln -sf ../../reads/hic_R1_001.fastq${GZ} omnic_R1_001.fastq${GZ}
ln -sf ../../reads/hic_R2_001.fastq${GZ} omnic_R2_001.fastq${GZ}
sbatch run_all.sh
cd ..

cd ..

###################################################################################################
## `01-asm`: contig assembly
###################################################################################################
cd 01-asm

cd hifiasm
HIFIASM=$(sbatch run_hifiasm.sh | cut -f 4 -d' ')
cd ..

cp -r template-purge-dups hifiasm-pd
cd hifiasm-pd
ln -sf ../hifiasm/hifi.hifiasm.bp.p_ctg.fasta contigs.fasta
PD_PLOT=$(sbatch -d afterany:${HIFIASM} run_purge_dups_plot.sh | cut -f 4 -d' ')
## NOTE: confirm the output plot and edit `L`, `M`, `U` in `run_purge_dups.sh` if necessary
PD=$(sbatch -d afterany:${PD_PLOT} run_purge_dups.sh | cut -f 4 -d' ')
cd ..

cd ..
