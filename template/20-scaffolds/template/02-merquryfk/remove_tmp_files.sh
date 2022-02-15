#!/bin/bash
shopt -s expand_aliases && source ~/.bashrc && set -e || exit 1
source ../../../config.sh

ml ${_FASTK}

Fastrm hifi.fastk scaffolds.fastk scaffolds.fastk.hifi.fastk
rm -rf tmp/
