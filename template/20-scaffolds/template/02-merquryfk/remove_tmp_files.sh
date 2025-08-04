#!/bin/bash
source ../../../config.sh
set -eu
ml ${_FASTK}
set -x

Fastrm hifi.fastk scaffolds.fastk scaffolds.fastk.hifi.fastk
rm -r $TMPDIR
