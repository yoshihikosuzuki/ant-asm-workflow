#!/bin/bash
source ../../config/aux.sh
eval $(parse_yaml ../../config/workflow.yaml)
eval ${shell_prefix}
set -eu
eval ${activate_fastk}
set -x

Fastrm hifi.fastk
rm -r tmp
