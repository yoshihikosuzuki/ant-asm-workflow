#!/bin/bash
#SBATCH -J omnic-partition
#SBATCH -o omnic-partition.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 128
#SBATCH --mem=500G
#SBATCH -t 24:00:00
shopt -s expand_aliases && source ~/.bashrc || exit 1
source ../../config.sh

CONTIGS=contigs.fasta
HAP1=contigs.hap1.fasta
HAP2=contigs.hap2.fasta
READS_1=omnic_R1_001.fastq
READS_2=omnic_R2_001.fastq
N_THREADS=128

_CONTIGS=$(basename ${CONTIGS} .gz)
_CONTIGS=${_CONTIGS%.*}

ml samtools bwa Other/arima_pipeline Other/seqkit

#samtools faidx ${CONTIGS}
#bwa index ${CONTIGS}

# Read mapping
for READS in ${READS_1} ${READS_2}; do
    _READS=$(basename ${READS} .gz)
    _READS=${_READS%.*}
    _OUT_PREFIX=${_CONTIGS}.${_READS}
    _OUT_BAM=${_OUT_PREFIX}.filtered.sorted.bam
    bwa mem -t${N_THREADS} -B8 ${CONTIGS} ${READS} |
        filter_five_end.pl |
        samtools view -@${N_THREADS} -b - |
        samtools sort -@${N_THREADS} -o ${_OUT_BAM}
    samtools index -@${N_THREADS} ${_OUT_BAM}

    # Extract reads for each haplotype
    for HAP in ${HAP1} ${HAP2}; do
        _HAP=$(basename ${HAP} .gz)
        _HAP=${_HAP%.*}
        _HAP=${_HAP#contigs.}
        _OUT_RNAMES=${_OUT_PREFIX}.${_HAP}.rnames
        _OUT_READS=${_READS}.${_HAP}.fastq
        samtools view -@${N_THREADS} ${_OUT_BAM} $(seqkit fx2tab -n ${HAP} | tr '\r\n' ' ') |
            cut -f1 > ${_OUT_RNAMES}
        seqkit -j ${N_THREADS} grep -f ${_OUT_RNAMES} ${READS} >${_OUT_READS}
    done
done


