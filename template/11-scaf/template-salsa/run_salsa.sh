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

CONTIGS=contigs.fasta
READS_1=omnic_R1_001.fastq
READS_2=omnic_R2_001.fastq
MIN_MAPQ=10
N_ITERATION=10
N_THREADS=128

CONTIGS_BWA_PREFIX=${CONTIGS}
OUT_PREFIX=${CONTIGS%.*}.${READS_1%%_R*}
SORTED_BAM=${OUT_PREFIX}.sorted.bam
OUT_BAM=${OUT_PREFIX}.dedup.sorted.bam
OUT_BED=${OUT_BAM/.bam/.bed}
OUT_SALSA=${OUT_PREFIX}_salsa

ml samtools bwa picard Other/SALSA/2.3

#samtools faidx ${CONTIGS}
#bwa index ${CONTIGS}

### Read mapping

# 1. VGP method
ml Other/arima_pipeline

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

# 2. Phase Genomics method
#bwa mem -t${N_THREADS} -5SP -B8 ${CONTIGS_BWA_PREFIX} ${READS_1} ${READS_2} |
#    samtools view -@${N_THREADS} -b -q ${MIN_MAPQ} -F 2316 - |
#    samtools sort -@${N_THREADS} -o ${SORTED_BAM}

### Deduplication
java -jar -Xmx500G -Djava.io.tmpdir=tmp/ $PICARD MarkDuplicates REMOVE_DUPLICATES=true I=${SORTED_BAM} O=${OUT_BAM} M=${OUT_BAM}.metrics ASSUME_SORT_ORDER=coordinate

### Scaffolding
bedtools bamtobed -i ${OUT_BAM} |
    sort -k 4 >${OUT_BED}
run_pipeline.py -a ${CONTIGS} -l ${CONTIGS}.fai -b ${OUT_BED} -e DNASE -o ${OUT_SALSA} -m yes -p yes -i ${N_ITERATION}

### Generate .assembly and .hic
ml Other/3d-dna

# 1. Use SALSA's script (that cannot be editable by Juicebox)
# 3d-dna-fasta2assembly ${OUT_SALSA}/scaffolds_FINAL.fasta >${OUT_SALSA}/scaffolds_FINAL.assembly
# convert.sh ${OUT_SALSA}
# cd ${OUT_SALSA} &&
#     ln -sf chromosome_sizes.tsv scaffolds_FINAL.chrom_sizes;
#     ln -sf salsa_scaffolds.hic scaffolds_FINAL.hic;
#     cd ..

# 2. Use 3D-DNA (that requires remapping of Omni-C reads but can be editable by Juicebox)
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
bwa index ${SCAFS}

3d-dna-fasta2assembly ${SCAFS} >${SCAFS%.*}.assembly
awk 'NF == 3 {print substr($1,2) "\t" $3}' ${SCAFS%.*}.assembly >${SCAFS%.*}.chrom_sizes

./scripts/juicer.sh -t ${N_THREADS} -S early -z ${SCAFS} -p ${SCAFS%.*}.chrom_sizes -g ${OUT_PREFIX}
mkdir -p hic &&
    cd hic &&
    3d-dna-run-assembly-visualizer ../${SCAFS%.*}.assembly ../aligned/merged_nodups.txt &&
    cd ..

ln -sf ${SCAFS%.*}.assembly .
ln -sf ${SCAFS%.*}.chrom_sizes .
ln -sf hic/scaffolds_FINAL.hic .
