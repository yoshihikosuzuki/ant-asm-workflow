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
