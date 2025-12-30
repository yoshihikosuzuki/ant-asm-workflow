#!/bin/bash
#SBATCH -J contig-eval
#SBATCH -o run_all.log
#SBATCH -p compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=1G
#SBATCH -t 72:00:00
set -euo pipefail

submit() { sbatch --parsable "$@"; }

MAKE_INDEX=$(submit 00-make_index.sh)

COMPLEASM=$( (cd 01-compleasm && submit run_compleasm.sh) )
MERQURY=$(   (cd 02-merqury   && submit run_merqury.sh) )

WINNOWMAP=$( (cd 04-winnowmap && submit --dependency=afterok:${MAKE_INDEX} run_winnowmap.sh) )
DEEPVAR=$(   (cd 05-deepvariant && submit --dependency=afterok:${WINNOWMAP} run_deepvariant.sh) )
MAPQV=$(     (cd 06-mapqv     && submit --dependency=afterok:${DEEPVAR} run_mapqv.sh) )

ASSET=$(     (cd 07-asset    && submit --dependency=afterok:${WINNOWMAP} run_asset.sh) )
TELOMERE=$(  (cd 09-telomere && submit run_make_telomere_bed.sh) )

# If you don't want Hi-C for contig eval, comment the next 2 lines out.
HIC=$(       (cd 11-hic      && submit --dependency=afterok:${MAKE_INDEX} run_hic.sh) )

FINAL_DEP="afterok:${COMPLEASM}:${MERQURY}:${MAPQV}:${ASSET}:${TELOMERE}:${HIC}"
FINAL=$(sbatch --parsable --dependency=${FINAL_DEP} -J contig-eval-done -o done.log --wrap "date; echo DONE")

echo "Submitted contig-eval. Final sentinel job: ${FINAL}"
