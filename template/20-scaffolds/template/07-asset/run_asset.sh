#!/bin/bash
#SBATCH -J asset
#SBATCH -o asset.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 128
#SBATCH --mem=100G
#SBATCH -t 24:00:00
shopt -s expand_aliases && source ~/.bashrc && set -e || exit 1
source ../../../config.sh

SCAF=scaffolds.fasta
READS=hifi.fastq
HIC_READS_1=omnic_R1_001.fastq
HIC_READS_2=omnic_R2_001.fastq
PB_BAM=scaffolds.hifi.winnowmap.sorted.bam
N_THREADS=128

OUT_PREFIX=${SCAF%.*}
SCAF_SPLIT=${OUT_PREFIX}.split.fasta
BED_SCAF=${OUT_PREFIX}.bed
BED_SCAF_GAP=${OUT_PREFIX}.gaps.bed
BED_SCAF_ACC=${OUT_PREFIX}.acc.bed
BED_SCAF_OK=${OUT_PREFIX}.reliable.bed
BED_SCAF_NG=${OUT_PREFIX}.unreliable.bed
PB_SAM=${PB_BAM/.bam/.sam}
PB_PAF=${PB_BAM/.bam/.paf}
BED_PB=${OUT_PREFIX}.pb.bed
BED_PB_OK=${OUT_PREFIX}.pb.reliable.bed
HIC_BAM=${OUT_PREFIX}.hic.bam
HIC_LINKS=${OUT_PREFIX}.hic.links.mat
BED_HIC=${OUT_PREFIX}.hic.bed
BED_HIC_OK=${OUT_PREFIX}.hic.reliable.bed

ml ${_ASSET}

# Contiguous regions and gaps
#samtools faix ${SCAF}
awk '{print $1"\t0\t"$2}' ${SCAF}.fai >${BED_SCAF}
detgaps ${SCAF} >${BED_SCAF_GAP}

bed_to_support() {
    IN_BED=$1
    _PREFIX=${IN_BED%.bed}
    BED_NG=${_PREFIX}.unreliable.bed
    BED_OK=${_PREFIX}.reliable.bed
    bedtools subtract -a ${BED_SCAF} -b ${IN_BED} |
        bedtools merge -d 100 -i - >${BED_NG}
    bedtools subtract -a ${BED_SCAF} -b ${BED_NG} >${BED_OK}
}

# PacBio
samtools view -@${N_THREADS} -h ${PB_BAM} -o ${PB_SAM}
paftools.js sam2paf ${PB_SAM} >${PB_PAF}
ast_pb -m3 ${PB_PAF} >${BED_PB}
bed_to_support ${BED_PB}

# Hi-C
split_fa ${SCAF} >${SCAF_SPLIT}
samtools faidx ${SCAF_SPLIT}
bwa index ${SCAF_SPLIT}
bwa mem -t${N_THREADS} -5SP -B8 ${SCAF_SPLIT} ${HIC_READS_1} ${HIC_READS_2} |
    samtools view -@${N_THREADS} -b -o ${HIC_BAM}
col_conts ${HIC_BAM} >${HIC_LINKS}
ast_hic2 ${SCAF_SPLIT}.fai ${HIC_LINKS} >${BED_HIC}
bed_to_support ${BED_HIC}

# Reliable block
acc ${BED_SCAF_GAP} ${BED_PB_OK} ${BED_HIC_OK} |
    awk '$4>=2' |
    bedtools merge -i - >${BED_SCAF_ACC}
bedtools subtract -a ${BED_SCAF} -b ${BED_SCAF_ACC} |
    bedtools merge -d 100 -i - >${BED_SCAF_NG}
bedtools subtract -a ${BED_SCAF} -b ${BED_SCAF_NG} >${BED_SCAF_OK}

touch reliable_blocks.n50
echo -n "Reliable block N50 length = " >reliable_blocks.n50
awk '{print $3 - $2}' ${BED_SCAF_OK} |
    sort -nr |
    awk '{ sum += $0; print $0, sum }' |
    tac |
    awk 'NR==1 { halftot=$2/2 } lastsize>halftot && $2<halftot { print $1 } { lastsize=$2 }' >>reliable_blocks.n50

# Generate .bed files for JBAT

bed_to_jbat() {
    CHROM_SIZES=$1
    BED_FILE=$2
    awk 'BEGIN {l=0} FNR == NR {offset[$1]=l; l+=$2; next} {printf "assembly\t" offset[$1]+$2 "\t" offset[$1]+$3; for(i=4;i<=NF;i++) printf "\t" $i; print ""}' ${CHROM_SIZES} ${BED_FILE}
}

bed_to_jbat ../${SCAF/.fasta/.chrom_sizes} ${BED_SCAF_GAP} > ${BED_SCAF_GAP/.bed/.JBAT.bed}
bed_to_jbat ../${SCAF/.fasta/.chrom_sizes} ${BED_SCAF_NG} > ${BED_SCAF_NG/.bed/.JBAT.bed}

if [ "$AUTO_DEL" = "true" ]; then
    source ./remove_tmp_files.sh
fi
