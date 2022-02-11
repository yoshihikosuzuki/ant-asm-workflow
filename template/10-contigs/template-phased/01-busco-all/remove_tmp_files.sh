#!/bin/bash
shopt -s expand_aliases && source ~/.bashrc && set -e || exit 1

rm -rf contigs.busco/*/ busco_downloads/ augustus_config/
