#!/bin/bash
#SBATCH -J omnic-partition
#SBATCH -o omnic-partition.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 128
#SBATCH --mem=500G
#SBATCH -t 24:00:00
source ../../config.sh
set -eu
module load ${_SAMTOOLS} ${_BWA} ${_ARIMA_PIPELINE} ${_SEQKIT} ${_PIGZ}
set -x

CONTIGS=contigs.fasta
HAP1=contigs.hap1.fasta
HAP2=contigs.hap2.fasta
READS_1=omnic_R1_001.fastq${OMNIC_GZ}
READS_2=omnic_R2_001.fastq${OMNIC_GZ}
N_THREADS=128

_CONTIGS=$(basename ${CONTIGS} .gz | sed 's/\.[^.]*$//')
READS_PREFIX=${READS_1%%_R*}

#samtools faidx ${CONTIGS}
#bwa index ${CONTIGS}

# Read mapping
for READS in ${READS_1} ${READS_2}; do
    _READS=$(basename ${READS} .gz | sed 's/\.[^.]*$//')
    _OUT_PREFIX=${_CONTIGS}.${_READS}
    _OUT_BAM=${_OUT_PREFIX}.filtered.sorted.bam
    bwa mem -t${N_THREADS} -B8 ${CONTIGS} ${READS} |
        filter_five_end.pl |
        samtools view -@${N_THREADS} -b - |
        samtools sort -@${N_THREADS} -o ${_OUT_BAM}
    samtools index -@${N_THREADS} ${_OUT_BAM}
done

# Extract reads for each haplotype
for HAP in ${HAP1} ${HAP2}; do
    _HAP=$(basename ${HAP} .gz | sed 's/\.[^.]*$//')
    _HAP=${_HAP#contigs.}
    for READS in ${READS_1} ${READS_2}; do
        _READS=$(basename ${READS} .gz | sed 's/\.[^.]*$//')
        _OUT_PREFIX=${_CONTIGS}.${_READS}
        _OUT_BAM=${_OUT_PREFIX}.filtered.sorted.bam
        _OUT_RNAMES=${_OUT_PREFIX}.${_HAP}.rnames
        samtools view -@${N_THREADS} ${_OUT_BAM} $(seqkit fx2tab -n ${HAP} | tr '\r\n' ' ') |
            cut -f1 |
            sort >${_OUT_RNAMES}
    done
    # Take union
    _HAP_RNAMES=${READS_PREFIX}.${_HAP}.rnames
    sort *.${_HAP}.rnames | uniq >${_HAP_RNAMES}

    for READS in ${READS_1} ${READS_2}; do
        _READS=$(basename ${READS} .gz | sed 's/\.[^.]*$//')
        _OUT_READS=${_READS}.${_HAP}.fastq.gz
        seqkit -j ${N_THREADS} grep -f ${_HAP_RNAMES} ${READS} | pigz -c >${_OUT_READS}
    done
done

if [ "$AUTO_DEL" = "true" ]; then
    source ./remove_tmp_files.sh
fi
