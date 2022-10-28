#!/bin/bash
source ../../config/aux.sh
eval $(parse_yaml ../../config/workflow.yaml)
eval ${shell_prefix}
set -eu
eval ${activate_genescope}
eval ${activate_merquryfk}
set -x

FASTK_PREFIX=hic.fastk
K=${genescope_hic_k}
PLOIDY=${genescope_ploidy}
HIST_MAX=1000
THRES_ERROR=${genescope_hic_thres_error}
N_THREAD=${genescope_threads}

OUT_PREFIX=${FASTK_PREFIX/.fastk/.genescope}

Histex -G ${FASTK_PREFIX} -h${HIST_MAX} |
    GeneScopeFK.R -o ${OUT_PREFIX} -p${PLOIDY} -k ${K}
KatGC -pdf -T${N_THREAD} -o${FASTK_PREFIX}.katgc ${FASTK_PREFIX}
PloidyPlot -v -pdf -T${N_THREAD} -e${THRES_ERROR} -o${FASTK_PREFIX}.ploidyplot ${FASTK_PREFIX}
