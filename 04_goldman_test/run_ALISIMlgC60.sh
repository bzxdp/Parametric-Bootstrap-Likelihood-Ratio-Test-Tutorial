#!/bin/bash
#SBATCH --job-name=ALISIM
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=100G
#SBATCH --time=24:00:00

cd $SLURM_SUBMIT_DIR

source ~/miniforge3/etc/profile.d/conda.sh
conda activate /projects/b56o/conda/envs/phylodev

# Alignment file: Nematoda.phy (see ../data/Nematoda.phy in this repo).
iqtree3 --alisim nematoda_sim \
    -s Nematoda.phy \
    -m LG+C60+F+G4{0.5408} \
    -te Nematoda.phy.treefile \
    --site-freq SAMPLING \
    --site-rate SAMPLING \
    --length 35371 \
    --num-alignments 100 \
    -af fasta \
    -seed 42 \
    -T 32 \
    --prefix nematedes_alisim 
