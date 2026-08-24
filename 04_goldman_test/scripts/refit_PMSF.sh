#!/bin/bash
# Re-fit each parametric-bootstrap replicate under CAT-PMSF, using the
# site-frequency profile fitted on the REAL data (03_pmsf_fit/output/nema2.sitefreq)
# frozen -- i.e. the profiles themselves are not re-estimated, only branch
# lengths / topology / gamma shape for that replicate.
#
# Usage: ./refit_PMSF.sh <sim_dir> <sitefreq_file> <out_dir> [n_threads]
set -euo pipefail
SIM_DIR=${1:?path to directory with nematoda_sim_*.fa}
SITEFREQ=${2:?path to nema2.sitefreq}
OUT_DIR=${3:?output directory}
NT=${4:-8}

mkdir -p "$OUT_DIR"
for fa in "$SIM_DIR"/nematoda_sim_*.fa; do
    base=$(basename "$fa" .fa)
    iqtree3 -s "$fa" -m Poisson+G4 -fs "$SITEFREQ" -nt "$NT" -pre "$OUT_DIR/${base}_pmsf_"
done
