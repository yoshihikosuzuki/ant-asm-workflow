#!/bin/bash
#SBATCH -J yahs
#SBATCH -o yahs.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 128
#SBATCH --mem=500G
#SBATCH -t 48:00:00
source ../../config.sh
set -eu
module load ${_SAMTOOLS} ${_BWA} ${_PICARD} ${_ARIMA_PIPELINE} ${_3DDNA} ${_YAHS} ${_SEQKIT}
set -x

CONTIGS=contigs.fasta
READS_1=omnic_R1_001.fastq${OMNIC_GZ}
READS_2=omnic_R2_001.fastq${OMNIC_GZ}

MIN_MAPQ=${SCAF_MIN_MAPQ}
ENZYME_NAME=${HIC_ENZYME_NAME}
N_THREADS=128

CONTIGS_BWA_PREFIX=${CONTIGS}
OUT_PREFIX=${CONTIGS%.*}.${READS_1%%_R*}
SORTED_BAM=${OUT_PREFIX}.sorted.bam
OUT_BAM=${OUT_PREFIX}.dedup.sorted.bam
OUT_SCAF=yahs.out_scaffolds_final.fa

# Read mapping
for READS in ${READS_1} ${READS_2}; do
    _READS=$(basename ${READS} .gz | sed 's/\.[^.]*$//')
    _OUT_BAM=${CONTIGS%.*}.${_READS}.filtered.bam
    bwa mem -t${N_THREADS} -B8 ${CONTIGS_BWA_PREFIX} ${READS} |
        filter_five_end.pl |
        samtools view -@${N_THREADS} -b -o ${_OUT_BAM}
done
two_read_bam_combiner.pl *.filtered.bam $(which samtools) ${MIN_MAPQ} |
    samtools view -@${N_THREADS} -b - |
    samtools sort -@${N_THREADS} -o ${SORTED_BAM}

# Deduplication
java -jar -Xmx500G -Djava.io.tmpdir=${TMPDIR} ${PICARD} MarkDuplicates REMOVE_DUPLICATES=true I=${SORTED_BAM} O=${OUT_BAM} M=${OUT_BAM}.metrics ASSUME_SORT_ORDER=coordinate

# YaHS scaffolding
mkdir output && mkdir output_for_curation
yahs ${CONTIGS} ${OUT_BAM}
cd output/
ln -sf ../${OUT_SCAF} scaffolds.fasta
seqkit stats -a scaffolds.fasta >scaffolds.fasta.stats
echo "Scaffold stats (output/scaffolds.fasta.stats):"
cat scaffolds.fasta.stats
cd ..

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

cd output_for_curation/
ln -sf ../references/scaffolds.fasta .
ln -sf ../references/scaffolds.assembly .
ln -sf ../hic/scaffolds.hic .
ln -sf ../aligned/merged_nodups.txt .
cd ..

if [ "$AUTO_DEL" = "true" ]; then
    source ./remove_tmp_files.sh
fi
