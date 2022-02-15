#!/bin/bash
shopt -s expand_aliases && source ~/.bashrc && set -e || exit 1
source ../../../config.sh

ml ${_FASTK}

Fastrm hifi.fastk contigs.hap1.fastk contigs.hap2.fastk contigs.hap1.fastk.hifi.fastk contigs.hap2.fastk.hifi.fastk
rm -rf tmp/
