#!/bin/bash
source ../../../config.sh
set -eu
module load ${_FASTK}
set -x

Fastrm hifi.fastk contigs.fastk contigs.fastk.hifi.fastk
rm -r $TMPDIR
