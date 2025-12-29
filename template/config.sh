#!/bin/bash
# NOTE: This script is loaded at the beginnig of each task.

# ----------------------------------------------------------------------------------------------- #
# Environment modules/Lmod settings
# NOTE: Write any commands necessary to load the modules listed just below
# ----------------------------------------------------------------------------------------------- #

hostname
date
source /apps/free/lmod/lmod/init/bash
module purge
module load bioinfo-ugrp-modules DebianMed/12.0

## Specifying ./tmp as the temporary directory, since /scratch sometimes does not have enough
## space left on Deigo.
export TMPDIR=tmp
mkdir -p $TMPDIR

# ----------------------------------------------------------------------------------------------- #
# List of the names and versions of the dependent modules
# NOTE: Variable names (left-hand side) must not be changed.
#       Right-hand side depends on your environment.
# ----------------------------------------------------------------------------------------------- #

## General tools
_SEQKIT=Other/seqkit/2.10.0
_PIGZ=pigz/2.6
_SAMTOOLS=samtools/1.12
_BCFTOOLS=Other/bcftools/1.15.1
_BEDTOOLS=bedtools/v2.29.2
_PICARD=picard/2.7.0

## For mapping; could influence results of any downstream analysis using them.
_BWA=bwa/0.7.17-3
_MINIMAP2=Other/minimap2/2.30
_GFATOOLS=Other/gfatools/0.5
_WINNOWMAP=Other/winnowmap/2.03

## For contig assembly
_HIFIASM=Other/hifiasm/0.25.0
_PURGE_DUPS=Other/purge_dups/1.2.5

## For scaffolding
_ARIMA_PIPELINE=Other/arima_pipeline/2019.02.08
_YAHS=Other/yahs/1.2.2

## For annotation and other analyses
_MAKE_TELOMERE_BED=Other/make_telomere_bed/2021.05.20
_MOSDEPTH=Other/mosdepth/0.3.11
_MUMMER=Other/mummer/4.0.0rc1

## For evaluation; could influence the quality values.
_BUSCO=Other/BUSCO/5.8.2
_COMPLEASM=Other/compleasm/0.2.7
_FASTK=Other/FASTK/2021.09.29
_GENESCOPE=Other/genescope/2021.03.26
_MERQURY=Other/merqury/1.3
_MERQURYFK=Other/MerquryFK/2021.09.14
_DEEPVARIANT=Other/deepvariant/1.9.0
_ASSET=Other/asset/1.0.3
_3DDNA=Other/3d-dna/180922

# ----------------------------------------------------------------------------------------------- #
# Automatic deletion of intermediate/temporary files
# - AUTO_DEL: Must be "true" or "false". If "true", some (large) intermediate files are 
#             removed immediately after each task finishes to save the disc space.
#             Even if you specify "false", you can manually delete the same files by
#             running `remove_tmp_files.sh` put in the directory of each task.
# ----------------------------------------------------------------------------------------------- #

AUTO_DEL=true

# ----------------------------------------------------------------------------------------------- #
# Whether or not input read files are gzipped or not.
# - HIFI_GZ, OMNIC_GZ: Set "" (empty string) if the HiFi/Omni-C read files are plain fastq files.
#                      Set ".gz" if gzipped.
# ----------------------------------------------------------------------------------------------- #

HIFI_GZ=".gz"
# HIFI_GZ=""
OMNIC_GZ=".gz"
# OMNIC_GZ=""

# ----------------------------------------------------------------------------------------------- #
# GeneScope/PloidyPlot (`00-data/hifi/` and `00-data/omnic/`) settings
# - PLOIDY: Ploidy of the genome.
# - HIFI_K: Length of k-mers used for HiFi reads by GeneScope/PoidyPlot.
# - OMNIC_K: Length of k-mers used for Omni-C reads by GeneScope/PoidyPlot.
# - HIFI_THRES_ERROR: K-mer counts smaller than this value are deemed as erroneous and excluded
#                     from PloidyPlot for HiFi reads. This does not have to be very accurate but
#                     should hopefully be decided based on the k-mer count histogram by e.g.
#                     Genescope.
# - OMNIC_THRES_ERROR: The same threshold for Omni-C reads.
# ----------------------------------------------------------------------------------------------- #

PLOIDY=2
HIFI_K=40
OMNIC_K=21
HIFI_THRES_ERROR=4
OMNIC_THRES_ERROR=4

# ----------------------------------------------------------------------------------------------- #
# Scaffolding (`11-scaf/template-yahs`) settings
# - SCAF_MIN_MAPQ: Minimum required MAPQ value (of BAM) for Omni-C read mappings used for scaffolding.
# - HIC_ENZYME_NAME: Name of the restriction enzyme used in the Hi-C experiment.
#                    Leave this empty for Omni-C.
# ----------------------------------------------------------------------------------------------- #

SCAF_MIN_MAPQ=10
HIC_ENZYME_NAME=

# ----------------------------------------------------------------------------------------------- #
# BUSCO/compleasm (`10-contigs/` and `20-scaffolds/`) settings
# - BUSCO_DB: Name of the BUSCO database. Depends on the species.
# - BUSCO_DB_DIR: Directory where you downloaded the DB with `compleasm download ${BUSCO_DB}`
# - USE_AUGUSTUS: Only for BUSCO. Must be "true" or "false". If "true", use Augustus (instead of
#                 Metaeuk) for gene prediction.
# ----------------------------------------------------------------------------------------------- #

BUSCO_DB="hymenoptera_odb12"
BUSCO_DB_DIR="/flash/EconomoU/gaurav/genome_asm/mb_downloads"
USE_AUGUSTUS=true

# ----------------------------------------------------------------------------------------------- #
# Merqury(FK) (`[10-contigs|20-scaffolds]/02-merqury(fk)`) settings
# - MERQURY_K: Length of k-mers used for Merqury(FK).
#              K ~ 20 should work for most genomes. The script `best_k.sh` provided by Merqury 
#              (e.g. `$ best_k.sh <GENOME_SIZE>`) will tell the appropiate value.
# ----------------------------------------------------------------------------------------------- #

MERQURY_K=19

# ----------------------------------------------------------------------------------------------- #
# Mapping-based QV (`[10-contigs|20-scaffolds]/06-mapqv`) settings
# - MAPQ_MIN_QUAL: Minimum required QUAL value (of GFF) for variants to be regarded as assembly
#                  errors. A larger value results in a higher mapping QV.
# - MAPQ_MIN_DEPTH, MAPQ_MAX_DEPTH: Variants whose depth (DP column of GFF) is within this range
#                                   are regarded as assembly errors. Should be changed if the
#                                   HiFi reads have a too high sequencing coverage.
# ----------------------------------------------------------------------------------------------- #

MAPQ_MIN_QUAL=1
MAPQ_MIN_DEPTH=5
MAPQ_MAX_DEPTH=200

# ----------------------------------------------------------------------------------------------- #
# Telomere (`20-scaffolds/09-telomere`) settings
# - TELOMERE_MOTIF: Unit sequence of the telomeric tandem array. Depends on the species.
# - TELOMERE_MIN_NCOPY: Minimum required copy number of the unit sequence in a tandem repeat array
#                       to be reported as telomeric.
# ----------------------------------------------------------------------------------------------- #

TELOMERE_MOTIF="TTAGG"
TELOMERE_MIN_NCOPY=100
