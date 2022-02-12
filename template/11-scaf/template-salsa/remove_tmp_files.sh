#!/bin/bash
shopt -s expand_aliases && source ~/.bashrc && set -e || exit 1

rm *.bam *.bed
cd *_salsa
rm -rf *.bed splits/ aligned/*.sam aligned/dups.txt aligned/merged_sort.txt aligned/opt_dups.txt hic/temp.*
