#!/bin/bash
shopt -s expand_aliases && source ~/.bashrc && set -e || exit 1
source ../../../config.sh

ml ${_FASTK}

Fastrm hifi.fastk contigs.fastk contigs.fastk.hifi.fastk
rm -rf tmp/
