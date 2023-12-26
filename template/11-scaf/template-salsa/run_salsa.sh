#!/bin/bash
#SBATCH -J salsa
#SBATCH -o salsa.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 128
#SBATCH --mem=500G
#SBATCH -t 48:00:00
shopt -s expand_aliases && source ~/.bashrc || exit 1
source ../../config.sh

mkdir output && mkdir output_for_curation

CONTIGS=contigs.fasta
READS_1=omnic_R1_001.fastq
READS_2=omnic_R2_001.fastq
MIN_MAPQ=${SALSA_MIN_MAPQ}
N_ITERATION=${SALSA_N_ITERATION}
ENZYME_NAME=${HIC_ENZYME_NAME}
N_THREADS=128

CONTIGS_BWA_PREFIX=${CONTIGS}
OUT_PREFIX=${CONTIGS%.*}.${READS_1%%_R*}
SORTED_BAM=${OUT_PREFIX}.sorted.bam
OUT_BAM=${OUT_PREFIX}.dedup.sorted.bam
OUT_BED=${OUT_BAM/.bam/.bed}
OUT_SALSA=${OUT_PREFIX}_salsa

ml ${_SAMTOOLS} ${_BWA} ${_PICARD} ${_ARIMA_PIPELINE} ${_SALSA} ${_SEQKIT}

# Read mapping
for READS in ${READS_1} ${READS_2}; do
    _OUT_PREFIX=${READS%.gz}
    _OUT_BAM=${CONTIGS%.*}.${_OUT_PREFIX%.*}.filtered.bam
    bwa mem -t${N_THREADS} -B8 ${CONTIGS_BWA_PREFIX} ${READS} |
        filter_five_end.pl |
        samtools view -@${N_THREADS} -b -o ${_OUT_BAM}
done
two_read_bam_combiner.pl *.filtered.bam $(which samtools) ${MIN_MAPQ} |
    samtools view -@${N_THREADS} -b - |
    samtools sort -@${N_THREADS} -o ${SORTED_BAM}
# Deduplication
java -jar -Xmx500G -Djava.io.tmpdir=tmp/ $PICARD MarkDuplicates REMOVE_DUPLICATES=true I=${SORTED_BAM} O=${OUT_BAM} M=${OUT_BAM}.metrics ASSUME_SORT_ORDER=coordinate

# Scaffolding
bedtools bamtobed -i ${OUT_BAM} | sort -k 4 >${OUT_BED}
run_pipeline.py -a ${CONTIGS} -l ${CONTIGS}.fai -b ${OUT_BED} -e DNASE -o ${OUT_SALSA} -m yes -p yes -i ${N_ITERATION}
cd output/ && ln -sf ../${OUT_SALSA}/scaffolds_FINAL.fasta ./scaffolds.fasta && cd ..

echo "Scaffold stats (${OUT_SALSA}/scaffolds_FINAL.fasta):"
seqkit stats -a ${OUT_SALSA}/scaffolds_FINAL.fasta

# Generate .assembly and .hic
ml ${_3DDNA}

cd ${OUT_SALSA}

# Make ./scripts/
juicer_copy_scripts_dir
# Remove ./aligned/ that already exists
rm -rf aligned
# Make ./fastq/
mkdir -p fastq && cd fastq && ln -sf ../../${READS_1} ../../${READS_2} . && cd ..
# Make ./references/
mkdir -p references && cd references && ln -sf ../scaffolds_FINAL.fasta . && cd ..

SCAFS=references/scaffolds_FINAL.fasta
CHROM_SIZES=${SCAFS%.*}.chrom_sizes
ASSEMBLY=${SCAFS%.*}.assembly

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

ln -sf ${ASSEMBLY} .
ln -sf ${CHROM_SIZES} .
ln -sf hic/scaffolds_FINAL.hic .

# Generate .mcool file
# IN_HIC=scaffolds_FINAL.hic
# OUT_COOL=scaffolds_FINAL.cool

# ml ${_HIC2COOL}

# hic2cool convert ${IN_HIC} ${OUT_COOL} -p ${N_THREADS}

cd ..

cd output_for_curation/
ln -sf ../${OUT_SALSA}/scaffolds_FINAL.fasta ./scaffolds.fasta
ln -sf ../${OUT_SALSA}/scaffolds_FINAL.assembly ./scaffolds.assembly
ln -sf ../${OUT_SALSA}/scaffolds_FINAL.hic ./scaffolds.hic
ln -sf ../${OUT_SALSA}/aligned/merged_nodups.txt .
cd ..

if [ "$AUTO_DEL" = "true" ]; then
    source ./remove_tmp_files.sh
fi
