#!/bin/bash
shopt -s expand_aliases && source ~/.bashrc && set -e || exit 1

rm -rf scaffolds.busco/*/ busco_downloads/ augustus_config/
