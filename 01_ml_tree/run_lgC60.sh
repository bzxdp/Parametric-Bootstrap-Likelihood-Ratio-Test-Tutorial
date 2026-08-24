#!/bin/bash
#SBATCH --job-name=iqtree
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=24:00:00

cd $SLURM_SUBMIT_DIR

source ~/miniforge3/etc/profile.d/conda.sh
conda activate /projects/b56o/conda/envs/phylodev

iqtree3 -s nematode2.phy -m LG+F+G+C60 -nt 16
