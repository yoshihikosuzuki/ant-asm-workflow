#!/bin/bash
shopt -s expand_aliases && source ~/.bashrc && set -e || exit 1

rm -rf splits/ aligned/*.sam aligned/dups.txt aligned/merged_sort.txt aligned/opt_dups.txt
