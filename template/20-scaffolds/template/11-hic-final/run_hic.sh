#!/bin/bash
#SBATCH -J hic
#SBATCH -o hic.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 128
#SBATCH --mem=500G
#SBATCH -t 48:00:00
shopt -s expand_aliases && source ~/.bashrc || exit 1
source ../../../config.sh

SCAFS=scaffolds.fasta
READS_1=omnic_R1_001.fastq
READS_2=omnic_R2_001.fastq
ENZYME_NAME=${HIC_ENZYME_NAME}
N_THREADS=128

OUT_PREFIX=${SCAFS%.*}

# Generate .assembly and .hic
ml ${_SAMTOOLS} ${_BWA} ${_3DDNA}

# Make ./scripts/
juicer_copy_scripts_dir
# Remove ./aligned/ that already exists
rm -rf aligned
# Make ./fastq/
mkdir -p fastq && cd fastq && ln -sf ../${READS_1} ../${READS_2} . && cd ..
# Make ./references/
mkdir -p references && cd references && ln -sf ../${SCAFS%.*}.* . && cd ..

SCAFS=references/${SCAFS}
CHROM_SIZES=${SCAFS%.*}.chrom_sizes
# ASSEMBLY=${SCAFS%.*}.assembly

# 3d-dna-fasta2assembly ${SCAFS} >${ASSEMBLY}
# awk 'NF == 3 {print substr($1,2) "\t" $3}' ${ASSEMBLY} >${CHROM_SIZES}

# samtools faidx ${SCAFS}
# bwa index ${SCAFS}

if [ -z "${ENZYME_NAME}" ]; then
    JUICER_S_OPT=""
else
    mkdir -p restriction_sites && cd restriction_sites
    generate_site_positions.py ${ENZYME_NAME} ${OUT_PREFIX} ../${SCAFS}
    cd ..
    JUICER_S_OPT="-s ${ENZYME_NAME}"
fi
./scripts/juicer.sh -t ${N_THREADS} -z ${SCAFS} -p ${CHROM_SIZES} -g ${OUT_PREFIX} ${JUICER_S_OPT}

if [ "$AUTO_DEL" = "true" ]; then
    source ./remove_tmp_files.sh
fi
