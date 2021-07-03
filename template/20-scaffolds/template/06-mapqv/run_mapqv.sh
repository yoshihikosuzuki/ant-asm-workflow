#!/bin/bash
#SBATCH -J mapqv
#SBATCH -o mapqv.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=500G
#SBATCH -t 24:00:00
shopt -s expand_aliases && source ~/.bashrc && set -e || exit 1

REF=scaffolds.fasta
READS=hifi.fastq
MIN_QUAL=30
MIN_DEPTH=5
MAX_DEPTH=200

# NOTE: Assuming the specific directory structure for input BAM and VCF files
_REF=$(basename ${REF} .gz)
_READS=$(basename ${READS} .gz)
IN_BAM=../04-winnowmap/${_REF%.*}.${_READS%.*}.winnowmap.sorted.bam
IN_VCF=../05-deepvariant/${_REF%.*}.${_READS%.*}.deepvariant.vcf

OUT_PREFIX=${IN_VCF%.vcf}
OUT_NORM_VCF=${OUT_PREFIX}.norm.vcf
OUT_VCF=${OUT_PREFIX}.norm.filtered.vcf
OUT_SNV_VCF=${OUT_PREFIX}.norm.filtered.snv.vcf

ml bcftools bedtools samtools

# Count the total number of bases
samtools view -F 0x100 -u ${IN_BAM} |
    bedtools genomecov -ibam - -split >${IN_BAM}.genomecov
awk -v l=${MIN_DEPTH} -v h=${MAX_DEPTH} '{if ($1=="genome" && $2>l && $2<h) {numbp += $3}} END {print numbp}' ${IN_BAM}.genomecov >${IN_BAM}.numbp
NUM_BP=$(cat ${IN_BAM}.numbp)
info echo "Total num. bases in mappable regions = $NUM_BP"

# Filter variants
bcftools norm -f ${REF} ${IN_VCF} -Ov > ${OUT_NORM_VCF}
bcftools view -i "QUAL>${MIN_QUAL} && ${MIN_DEPTH}<DP && DP<${MAX_DEPTH} && (GT=\"AA\" || GT=\"Aa\")" -Ov ${OUT_NORM_VCF} >${OUT_VCF}
bcftools view -v snps -Ov ${OUT_VCF} >${OUT_SNV_VCF}

calc_stats() {
    VCF=$1
    # Count the total number of variants
    bcftools view -H -Ov ${VCF} |
        awk -F "\t" '{print $4"\t"$5}' |
        awk '{lenA=length($1); lenB=length($2); if (lenA < lenB) {sum+=lenB-lenA} else if (lenA > lenB) { sum+=lenA-lenB } else {sum+=lenA}} END {print sum}' >${VCF}.numvar
    NUM_VAR=$(cat ${VCF}.numvar)
    #info echo "Total num. bases subject to change: $NUM_VAR"
    # Calculate QV
    QV=$(echo "$NUM_VAR $NUM_BP" | awk '{print (-10*log($1/$2)/log(10))}')
    info echo -e "mapping QV (${VCF}) = $QV"
}

# Calculate QV
calc_stats ${OUT_VCF}
calc_stats ${OUT_SNV_VCF}
