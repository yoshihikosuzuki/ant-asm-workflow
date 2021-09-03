#!/bin/bash

### -------- Clean up and Load modules -------- ###

module purge
module use /apps/.bioinfo-ugrp-modulefiles81
module use /apps/unit/BioinfoUgrp/DebianMed/10.7/modulefiles


### -------- AUTOMATIC DELETION OF INTERMEDIATE FILES (every step) -------- ###

# NOTE: "true" or "false"
#       Even if you specify "false", you can manually delete the same files
#       by running `remove_tmp_files.sh` put in each directory.
AUTO_DEL=true


### -------- HiCanu (`01-asm/hicanu/`) -------- ###

# NOTE: This does not have to be very accurate.
GENOME_SIZE=300000000


### -------- 3D-DNA (`11-scaf/template-3ddna`) -------- ###

# NOTE: Leave this empty for Omni-C.
HIC_ENZYME_NAME=


### -------- BUSCO (`10-contigs/` and `20-scaffolds/`) -------- ###

BUSCO_DB="hymenoptera_odb10"


### -------- Merqury (`10-contigs/` and `20-scaffolds/`) -------- ###

# NOTE: K ~ 20 should work for most datasets, but
#       `$ ml Other/merqury; best_k.sh <GENOME_SIZE>` will tell the appropiate value.
MERQURY_K=19


### -------- Telomere (`20-scaffolds/`) -------- ###

TELOMERE_MOTIF="TTAGG"

