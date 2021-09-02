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
source ../../config.sh

CONTIGS=contigs.fasta
READS_1=omnic_R1_001.fastq
READS_2=omnic_R2_001.fastq
ENZYME_NAME=${HIC_ENZYME_NAME}
N_THREADS=128

OUT_PREFIX=${CONTIGS%.*}

ml Other/3d-dna Other/seqkit

# Make ./scripts/
juicer_copy_scripts_dir
# Remove ./aligned/ that already exists
rm -rf aligned
# Make ./fastq/
mkdir -p fastq && cd fastq && ln -sf ../${READS_1} ../${READS_2} . && cd ..
# Make ./references/
mkdir -p references && cd references && ln -sf ../${CONTIGS}* . && cd ..

CONTIGS=references/$(basename ${CONTIGS})
CHROM_SIZES=${CONTIGS%.*}.chrom_sizes
ASSEMBLY=${CONTIGS%.*}.assembly

#bwa index ${CONTIGS}

# Generate .assembly and .chrom_sizes files
3d-dna-fasta2assembly ${CONTIGS} >${ASSEMBLY}
awk 'NF == 3 {print substr($1,2) "\t" $3}' ${ASSEMBLY} >${CHROM_SIZES}

# Generate `merged_nodups.txt`
if [ -z "${ENZYME_NAME}" ]; then
    JUICER_S_OPT=""
else
    mkdir -p restriction_sites && cd restriction_sites
    generate_site_positions.py ${ENZYME_NAME} ${OUT_PREFIX} ../${CONTIGS}
    cd ..
    JUICER_S_OPT="-s ${ENZYME_NAME}"
fi
./scripts/juicer.sh -t ${N_THREADS} -S early -z ${CONTIGS} -p ${CHROM_SIZES} -g ${OUT_PREFIX} ${JUICER_S_OPT}

# Generate .hic file for contigs before scaffolding
# mkdir -p hic && cd hic
# 3d-dna-run-assembly-visualizer ../${ASSEMBLY} ../aligned/merged_nodups.txt
# cd ..

# Scaffolding
# NOTE: Use `*.final.assembly` (not `.FINAL.assembly`) for visualization with `*.final.hic`
rm -rf scaffolding && mkdir -p scaffolding && cd scaffolding
ln -sf ../${CONTIGS} ../aligned/merged_nodups.txt .
3d-dna ${CONTIGS##*/} merged_nodups.txt
cd ..

# Generate .mcool file
# NOTE: HiGlass is not editable, so `.FINAL.assembly` should be preferable for `.chrom_sizes`.

IN_ASSEMBLY=contigs.FINAL.assembly
OUT_CHROM_SIZES=contigs.final.chrom_sizes
IN_HIC=contigs.final.hic
OUT_COOL=contigs.final.cool

ml Other/hic2cool

cd scaffolding
hic2cool convert ${IN_HIC} ${OUT_COOL} -p ${N_THREADS}
awk 'NF == 3 {print substr($1,2) "\t" $3}' ${IN_ASSEMBLY} >${OUT_CHROM_SIZES}
cd ..

echo "Scaffold stats (scaffolding/contigs.FINAL.fasta):"
seqkit stats -a scaffolding/contigs.FINAL.fasta

if [ "$AUTO_DEL" = "true" ]; then
    source ./remove_tmp_files.sh
fi
