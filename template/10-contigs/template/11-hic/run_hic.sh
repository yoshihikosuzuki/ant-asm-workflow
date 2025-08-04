#!/bin/bash
#SBATCH -J hic
#SBATCH -o hic.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 128
#SBATCH --mem=500G
#SBATCH -t 48:00:00
source ../../../config.sh
set -eu
module load ${_SAMTOOLS} ${_BWA} ${_3DDNA}
set -x

SCAFS=contigs.fasta
READS_1=omnic_R1_001.fastq${OMNIC_GZ}
READS_2=omnic_R2_001.fastq${OMNIC_GZ}
ENZYME_NAME=${HIC_ENZYME_NAME}
N_THREADS=128

OUT_PREFIX=${SCAFS%.*}

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
ASSEMBLY=${SCAFS%.*}.assembly

# NOTE: Index files should be generated before
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
./scripts/juicer.sh -t ${N_THREADS} -S early -z ${SCAFS} -p ${CHROM_SIZES} -g ${OUT_PREFIX} ${JUICER_S_OPT}

mkdir -p hic
cd hic
mkdir -p ${TMPDIR}
3d-dna-run-assembly-visualizer ../${ASSEMBLY} ../aligned/merged_nodups.txt || true
cd ..

if [ "$AUTO_DEL" = "true" ]; then
    source ./remove_tmp_files.sh
fi
