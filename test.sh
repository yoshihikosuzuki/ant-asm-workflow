#!/bin/bash
set -eux

PATH_TO_TEMPLATE=./template


rm -rf test
cp -r ${PATH_TO_TEMPLATE} test
cd test

# download test data
wget -O - https://mlab.cb.k.u-tokyo.ac.jp/~yoshihiko_s/ant-asm-workflow/reads.tar.gz | tar xzvf -

# `00-data`: make symlink to reads of the test data
cd 00-data

cd hifi
ln -sf ../../reads/hifi.fastq.gz .
cd ..

cd omnic
ln -sf ../../reads/hic_R1_001.fastq.gz omnic_R1_001.fastq.gz
ln -sf ../../reads/hic_R2_001.fastq.gz omnic_R2_001.fastq.gz
cd ..
