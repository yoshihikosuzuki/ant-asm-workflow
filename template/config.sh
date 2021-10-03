#!/bin/bash
# NOTE: This script is loaded at the beginnig of each task.

### -------- Clean up and Load modules -------- ###

module purge
module use /apps/.bioinfo-ugrp-modulefiles81
module use /apps/unit/BioinfoUgrp/DebianMed/10.7/modulefiles


### -------- AUTOMATIC DELETION OF INTERMEDIATE FILES (every step) -------- ###

# NOTE: Must be "true" or "false". If "true", some (large) intermediate files are 
#       removed immediately after each task finishes to save the disc space.
#       Even if you specify "false", you can manually delete the same files by
#       running `remove_tmp_files.sh` put in the directory of each task.
AUTO_DEL=true


### -------- Genescope/PloidyPlot (`00-data/hifi/` and `00-data/omnic/`) -------- ###

PLOIDY=2
HIFI_K=40
OMNIC_K=21
# NOTE: K-mer counts **smaller** than the threshold below are deemed as erroneous and
#       excluded from PloidyPlot. This does not have to be very accurate but should
#       hopefully be determined based on the k-mer count histogram (i.e. Gen(om)escope).
HIFI_THRES_ERROR=4
OMNIC_THRES_ERROR=4


### -------- HiCanu (`01-asm/hicanu/`) -------- ###

# NOTE: This does not have to be very accurate.
GENOME_SIZE=300000000


### -------- SALSA (`11-scaf/template-salsa`) -------- ###

SALSA_MIN_MAPQ=10
SALSA_N_ITERATION=10


### -------- 3D-DNA (`11-scaf/template-3ddna`) -------- ###

# NOTE: Leave this empty for Omni-C.
HIC_ENZYME_NAME=


### -------- BUSCO (`10-contigs/` and `20-scaffolds/`) -------- ###

BUSCO_DB="hymenoptera_odb10"


### -------- Merqury (`10-contigs/` and `20-scaffolds/`) -------- ###

# NOTE: K ~ 20 should work for most genomes. The script `best_k.sh` provided by
#       Merqury (e.g. `$ best_k.sh <GENOME_SIZE>`) will tell the appropiate value.
MERQURY_K=19


### -------- Telomere (`20-scaffolds/`) -------- ###

TELOMERE_MOTIF="TTAGG"
# NOTE: Min copy numer of the motif sequence to be reported
TELOMERE_MIN_NCOPY=100


### -------- NAMES OF ENVIRONMENT MODULES -------- ###

_SEQKIT=Other/seqkit/2.0.0
_SAMTOOLS=samtools/1.12
_BCFTOOLS=bcftools/1.9-1
_BEDTOOLS=bedtools/v2.29.2
_GFATOOLS=Other/gfatools/0.5
_PICARD=picard/2.7.0
_BWA=bwa/0.7.17-3
_MINIMAP2=Other/minimap2/2.20
_WINNOWMAP=Other/winnowmap/2.03
_HIFIASM=Other/hifiasm/0.15.4
_CANU=Other/canu/2.1.1
_IPA=Other/pbipa/1.3.2
_PEREGRINE=Other/peregrine/1.6.3
_PURGE_DUPS=Other/purge_dups/1.2.5
_ARIMA_PIPELINE=Other/arima_pipeline/2019.02.08
_SALSA=Other/SALSA/2.3
_3DDNA=Other/3d-dna/180922
_HIC2COOL=Other/hic2cool/0.8.3
_BUSCO=Other/BUSCO/5.1.3
_FASTK=Other/FASTK/2021.05.27
_GENESCOPE=Other/genescope/2021.03.26
_MERQURYFK=Other/MerquryFK/2021.09.14
_DEEPVARIANT=Other/deepvariant/1.1.0
_ASSET=Other/asset/1.0.3
_MAKE_TELOMERE_BED=Other/make_telomere_bed/2021.05.20
_MOSDEPTH=Other/mosdepth/0.3.1
