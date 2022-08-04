#!/bin/bash
#SBATCH -J yahs
#SBATCH -o yahs.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 128
#SBATCH --mem=500G
#SBATCH -t 48:00:00
shopt -s expand_aliases && source ~/.bashrc || exit 1
source ../../config.sh
ml ${_SAMTOOLS} ${_BWA} ${_3DDNA} ${_YAHS} ${_SEQKIT}

CONTIGS=contigs.fasta
BAM=contigs.omnic.dedup.sorted.bam
READS_1=omnic_R1_001.fastq
READS_2=omnic_R2_001.fastq
ENZYME_NAME=${HIC_ENZYME_NAME}
OUT_SCAF=yahs.out_scaffolds_final.fa

yahs ${CONTIGS} ${BAM}
cd output/ && ln -sf ../${OUT_SCAF} scaffolds.fasta && cd ..

echo "Scaffold stats (output/scaffolds.fasta):"
seqkit stats -a output/scaffolds.fasta

# Generate .assembly and .hic
# Make ./scripts/
juicer_copy_scripts_dir
# Remove ./aligned/ that already exists
rm -rf aligned
# Make ./fastq/
mkdir -p fastq && cd fastq && ln -sf ../${READS_1} ../${READS_2} . && cd ..
# Make ./references/
mkdir -p references && cd references && ln -sf ../output/scaffolds.fasta . && cd ..

SCAFS=references/scaffolds.fasta
CHROM_SIZES=${SCAFS%.*}.chrom_sizes
ASSEMBLY=${SCAFS%.*}.assembly
OUT_PREFIX=scaffolds.omnic
N_THREADS=128

samtools faidx ${SCAFS}
bwa index ${SCAFS}

# Generate .assembly and .chrom_sizes files
3d-dna-fasta2assembly ${SCAFS} >${ASSEMBLY}
awk 'NF == 3 {print substr($1,2) "\t" $3}' ${ASSEMBLY} >${CHROM_SIZES}

# Generate `merged_nodups.txt`
if [ -z "${ENZYME_NAME}" ]; then
    JUICER_S_OPT=""
else
    mkdir -p restriction_sites && cd restriction_sites
    generate_site_positions.py ${ENZYME_NAME} ${OUT_PREFIX} ../${SCAFS}
    cd ..
    JUICER_S_OPT="-s ${ENZYME_NAME}"
fi
./scripts/juicer.sh -t ${N_THREADS} -S early -z ${SCAFS} -p ${CHROM_SIZES} -g ${OUT_PREFIX} ${JUICER_S_OPT}

# Generate .hic file
mkdir -p hic && cd hic
3d-dna-run-assembly-visualizer ../${ASSEMBLY} ../aligned/merged_nodups.txt
cd ..

# ln -sf ${ASSEMBLY} .
# ln -sf ${CHROM_SIZES} .
# ln -sf hic/scaffolds_FINAL.hic .

# Generate .mcool file
# IN_HIC=scaffolds_FINAL.hic
# OUT_COOL=scaffolds_FINAL.cool

# ml ${_HIC2COOL}

# hic2cool convert ${IN_HIC} ${OUT_COOL} -p ${N_THREADS}

cd output_for_curation/
ln -sf ../references/scaffolds.fasta .
ln -sf ../references/scaffolds.assembly .
ln -sf ../hic/scaffolds.hic .
ln -sf ../aligned/merged_nodups.txt .
cd ..

if [ "$AUTO_DEL" = "true" ]; then
    source ./remove_tmp_files.sh
fi
