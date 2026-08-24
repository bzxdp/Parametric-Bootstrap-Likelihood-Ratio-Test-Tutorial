#!/bin/bash
# Re-fit each parametric-bootstrap replicate under the null model,
# LG+F+G+C60, allowing all free parameters (topology, branch lengths,
# C60 component weights, gamma shape) to be re-estimated.
#
# Usage: ./refit_LGC60.sh <sim_dir> <out_dir> [n_threads]
set -euo pipefail
SIM_DIR=${1:?path to directory with nematoda_sim_*.fa}
OUT_DIR=${2:?output directory}
NT=${3:-8}

mkdir -p "$OUT_DIR"
for fa in "$SIM_DIR"/nematoda_sim_*.fa; do
    base=$(basename "$fa" .fa)
    iqtree3 -s "$fa" -m LG+F+G+C60 -nt "$NT" -pre "$OUT_DIR/${base}_tree"
done
