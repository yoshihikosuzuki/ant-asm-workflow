#!/bin/bash
shopt -s expand_aliases && source ~/.bashrc && set -e || exit 1
source ../../../config.sh

ml ${_FASTK}

Fastrm hifi.fastk scaffolds.hap1.fastk scaffolds.hap2.fastk scaffolds.hap1.fastk.hifi.fastk scaffolds.hap2.fastk.hifi.fastk
rm -rf tmp/
