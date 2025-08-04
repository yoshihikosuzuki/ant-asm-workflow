#!/bin/bash
#SBATCH -J deepvariant
#SBATCH -o deepvariant.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 128
#SBATCH --mem=500G
#SBATCH -t 48:00:00
source ../../../config.sh
set -eu
module load ${_DEEPVARIANT}
set -x

REF=contigs.fasta
READS=hifi.fastq${HIFI_GZ}
IN_SORTED_BAM=contigs.hifi.winnowmap.sorted.bam
N_THREADS=128

_READS=$(basename ${READS} .gz | sed 's/\.[^.]*$//')
OUT_PREFIX=${REF%.*}.${_READS%.*}.deepvariant
OUT_VCF=${OUT_PREFIX}.vcf

run_deepvariant --num_shards ${N_THREADS} --model_type PACBIO --ref ${REF} --reads ${IN_SORTED_BAM} --output_vcf ${OUT_VCF} --intermediate_results_dir ${TMPDIR}

if [ "$AUTO_DEL" = "true" ]; then
    source ./remove_tmp_files.sh
fi
