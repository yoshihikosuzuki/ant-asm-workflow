#!/bin/bash

ROOT=$(readlink -f $(dirname $0))

eval "snakemake \
    -s ${ROOT}/main.smk \
    -d ${ROOT} \
    --profile ${ROOT}/config/ \
    --configfiles ${ROOT}/config/workflow.yaml \
    $@"
