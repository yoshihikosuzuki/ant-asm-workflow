#!/bin/bash
source ../../config/aux.sh
eval $(parse_yaml ../../config/workflow.yaml)
eval ${shell_prefix}
set -eu
eval ${activate_fastk}
set -x

IN_HIC="${input_hic}"
K=${genescope_hic_k}
N_THREAD=${fastk_threads}
N_MEMORY=${fastk_mem_gb}
TMP_DIR=tmp

OUT_PREFIX=hic.fastk

mkdir -p ${TMP_DIR}
FastK -k${K} -T${N_THREAD} -M${N_MEMORY} -v -t1 -p -P${TMP_DIR} -N${OUT_PREFIX} ${IN_HIC}
Tabex ${OUT_PREFIX} CHECK
