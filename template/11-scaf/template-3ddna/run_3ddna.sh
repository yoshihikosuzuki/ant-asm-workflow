#!/bin/bash
#SBATCH -J 3ddna
#SBATCH -o 3ddna.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 128
#SBATCH --mem=500G
#SBATCH -t 48:00:00
shopt -s expand_aliases && source ~/.bashrc || exit 1

CONTIGS=contigs.fasta
READS_1=omnic_R1_001.fastq
READS_2=omnic_R2_001.fastq
N_THREADS=128

OUT_PREFIX=${CONTIGS%.*}.

ml bwa 3d-dna

# Make ./scripts/
juicer_copy_scripts_dir
# Remove ./aligned/ that already exists
rm -rf aligned
# Make ./fastq/
mkdir -p fastq && cd fastq && ln -sf ../${READS_1} ../${READS_2} . && cd ..
# Make ./references/
mkdir -p references && cd references && ln -sf ../${CONTIGS}* . && cd ..

CONTIGS=references/${CONTIGS##*/}
#bwa index ${CONTIGS}

# Generate .assembly and .chrom_sizes files
3d-dna-fasta2assembly ${CONTIGS} >${CONTIGS%.*}.assembly
awk 'NF == 3 {print substr($1,2) "\t" $3}' ${CONTIGS%.*}.assembly >${CONTIGS%.*}.chrom_sizes

# Generate .hic file
./scripts/juicer.sh -t ${N_THREADS} -S early -z ${CONTIGS} -p ${CONTIGS%.*}.chrom_sizes -g ${OUT_PREFIX}
mkdir -p hic && cd hic &&
    3d-dna-run-assembly-visualizer ../${CONTIGS%.*}.assembly ../aligned/merged_nodups.txt
cd ..

# Scaffolding
# NOTE: Use `*.final.assembly` (not `.FINAL.assembly`) for visualization with `*.final.hic`
rm -rf scaffolding && mkdir -p scaffolding && cd scaffolding &&
    ln -sf ../${CONTIGS} ../aligned/merged_nodups.txt . &&
    3d-dna ${CONTIGS##*/} merged_nodups.txt
cd ..
