#!/bin/bash
source ../../config/aux.sh
eval $(parse_yaml ../../config/workflow.yaml)
eval ${shell_prefix}
set -eu
eval ${activate_fastk}
set -x

IN_HIFI=${input_hifi}
K=${genescope_hifi_k}
N_THREAD=${fastk_threads}
MEM_GB=${fastk_mem_gb}
TMP_DIR=tmp

_IN_HIFI=$(basename ${IN_HIFI} .gz)
OUT_PREFIX=${_IN_HIFI%.*}.fastk

mkdir -p ${TMP_DIR}
FastK -k${K} -T${N_THREAD} -M${MEM_GB} -v -t1 -p -P${TMP_DIR} -N${OUT_PREFIX} ${IN_HIFI}
Tabex ${OUT_PREFIX} CHECK
