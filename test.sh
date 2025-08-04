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

## Wait for the completion of the jobs
srun -p compute -c 1 --mem 1G -t 1:00:00 -d afterany:${PD} --wait=0 sleep 1s

###################################################################################################
## `10-contigs`: contig assembly evaluation
###################################################################################################
cd 10-contigs

cp -r template hifiasm-pd
cd hifiasm-pd
ln -sf ../../01-asm/hifiasm-pd/purged.fa contigs.fasta
bash run_all.sh
cd ..

cd ..

###################################################################################################
## `11-scaf`: scaffold assembly
###################################################################################################
cd 11-scaf

cp -r template-yahs hifiasm-pd-yahs
cd hifiasm-pd-yahs
ln -sf ../../10-contigs/hifiasm-pd/contigs.fasta* .
YAHS=$(sbatch run_yahs.sh | cut -f 4 -d' ')
cd ..

cd ..

## Wait for the completion of the jobs
srun -p compute -c 1 --mem 1G -t 1:00:00 -d afterany:${YAHS} --wait=0 sleep 1s

###################################################################################################
## `20-scaffolds`: scaffold assembly evaluation
###################################################################################################
cd 20-scaffolds

cp -r template hifiasm-pd-yahs
cd hifiasm-pd-yahs
ln -sf ../../11-scaf/hifiasm-pd-yahs/output/scaffolds.fasta scaffolds.fasta
bash run_all.sh
cd ..

cd ..

###################################################################################################

cd ..   # test
