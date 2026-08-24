#!/bin/bash

#SBATCH --job-name=pmsf_O
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=60GB
#SBATCH --time=24:00:00

source ~/miniforge3/etc/profile.d/conda.sh
conda activate /projects/b56o/conda/envs/phylodev

cd $SLURM_SUBMIT_DIR

# Alignment file: Nematoda.phy (see ../data/Nematoda.phy in this repo).
iqtree  -s Nematoda.phy  -m Poisson+G4 -fs nema2.sitefreq   -nt 32 -pre CAT_pmsf_
