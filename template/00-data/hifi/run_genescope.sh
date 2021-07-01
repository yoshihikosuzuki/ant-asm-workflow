#!/bin/bash
#SBATCH -J FASTK
#SBATCH -o fastk.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=500G
#SBATCH -t 24:00:00

K=
N_THREAD=16
N_MEMORY=16

ml FASTK
mkdir -p tmp

# Single dataset

IN_FNAMES=(   # <source_fname> per line

)

for IN_FNAME in "${IN_FNAMES[@]}"; do
    OUT_PREFIX=${IN_FNAME%.*}.fastk
    FastK -k${K} -T${N_THREAD} -M${N_MEMORY} -v -t1 -p -Ptmp -N${OUT_PREFIX} ${IN_FNAME}
    Tabex ${OUT_PREFIX} CHECK
done

# Relative profile

IN_FNAMES=(   # "<source_fname> <table_prefix>" per line

)

for INPUTS in "${IN_FNAMES[@]}"; do
    ARRAY=($INPUTS)
    IN_FNAME=${ARRAY[0]}
    IN_TABLE_PREFIX=${ARRAY[1]}
    FastK -k${K} -T${N_THREAD} -M${N_MEMORY} -v -p:${IN_TABLE_PREFIX} -Ptmp -N${IN_FNAME%.*}.${IN_TABLE_PREFIX} ${IN_FNAME}
done
