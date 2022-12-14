#!/bin/bash
# NOTE: This script is loaded at the beginnig of each task.

# ----------------------------------------------------------------------------------------------- #
# Environment modules/Lmod settings
# NOTE: Write any commands necessary to load the modules listed just below
# ----------------------------------------------------------------------------------------------- #

module purge
module use --append /apps/.modulefiles72
module use /apps/.bioinfo-ugrp-modulefiles81
module use /apps/unit/BioinfoUgrp/DebianMed/10.7/modulefiles

# ----------------------------------------------------------------------------------------------- #
# List of the names and versions of the dependent modules
# NOTE: Variable names (left-hand side) must not be changed.
#       Right-hand side depends on your environment.
# ----------------------------------------------------------------------------------------------- #

_SEQKIT=Other/seqkit/2.0.0
_SAMTOOLS=samtools/1.12
_BCFTOOLS=bcftools/1.9-1
_BEDTOOLS=bedtools/v2.29.2
_GFATOOLS=Other/gfatools/0.5
_PICARD=picard/2.7.0
_BWA=bwa/0.7.17-3
_MINIMAP2=Other/minimap2/2.20
_WINNOWMAP=Other/winnowmap/2.03
_HIFIASM=Other/hifiasm/0.16.1
_CANU=Other/canu/2.1.1
_IPA=Other/pbipa/1.3.2
_PEREGRINE=Other/peregrine/1.6.3
_PURGE_DUPS=Other/purge_dups/1.2.5
_ARIMA_PIPELINE=Other/arima_pipeline/2019.02.08
_SALSA=Other/SALSA/2.3
_3DDNA=Other/3d-dna/180922
_YAHS=Other/yahs/1.1
_HIC2COOL=Other/hic2cool/0.8.3
_BUSCO=Other/BUSCO/5.1.3
_FASTK=Other/FASTK/2021.09.29
_GENESCOPE=Other/genescope/2021.03.26
_MERQURY=Other/merqury/1.3
_MERQURYFK=Other/MerquryFK/2021.09.14
_MERYL=Other/meryl/1.3
_DEEPVARIANT=Other/deepvariant/1.1.0
_ASSET=Other/asset/1.0.3
_MAKE_TELOMERE_BED=Other/make_telomere_bed/2021.05.20
_MOSDEPTH=Other/mosdepth/0.3.1
_MUMMER=Other/mummer/4.0.0rc1

# ----------------------------------------------------------------------------------------------- #
# Automatic deletion of intermediate/temporary files
# - AUTO_DEL: Must be "true" or "false". If "true", some (large) intermediate files are 
#             removed immediately after each task finishes to save the disc space.
#             Even if you specify "false", you can manually delete the same files by
#             running `remove_tmp_files.sh` put in the directory of each task.
# ----------------------------------------------------------------------------------------------- #

AUTO_DEL=true

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
# HiCanu (`01-asm/hicanu/`) settings
# - GENOME_SIZE: Genome size in bp. This does not have to be very accurate.
# ----------------------------------------------------------------------------------------------- #

GENOME_SIZE=300000000

# ----------------------------------------------------------------------------------------------- #
# SALSA (`11-scaf/template-salsa`) settings
# - SALSA_MIN_MAPQ: Minimum required MAPQ value (of BAM) for Omni-C read mappings used for SALSA.
# - SALSA_N_ITERATION: Maximum number of iterartions of SALSA procedure.
# ----------------------------------------------------------------------------------------------- #

SALSA_MIN_MAPQ=10
SALSA_N_ITERATION=10

# ----------------------------------------------------------------------------------------------- #
# 3D-DNA (`11-scaf/template-3ddna`) settings
# - HIC_ENZYME_NAME: Name of the restriction enzyme used in the Hi-C experiment.
#                    Leave this empty for Omni-C.
# ----------------------------------------------------------------------------------------------- #

HIC_ENZYME_NAME=

# ----------------------------------------------------------------------------------------------- #
# BUSCO (`10-contigs/` and `20-scaffolds/`) settings
# - BUSCO_DB: Name of the BUSCO database. Depends on the species.
# - USE_AUGUSTUS: Must be "true" or "false". If "true", use Augustus (instead of Metaeuk) for
#                 gene prediction.
# ----------------------------------------------------------------------------------------------- #

BUSCO_DB="hymenoptera_odb10"
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
