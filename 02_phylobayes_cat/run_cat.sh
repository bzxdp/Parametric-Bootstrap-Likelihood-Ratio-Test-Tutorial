#!/bin/bash

#SBATCH --job-name=N1C
#SBATCH --nodes=3
#SBATCH --ntasks=96
#SBATCH --mem=60GB
#SBATCH --time=24:00:00

module load gnu12/12.2.0
module load openmpi4/4.1.5

export OMPI_MCA_pml=ob1
export OMPI_MCA_btl=self,tcp


cd $SLURM_SUBMIT_DIR


srun -n 96 /projects/b56o/software/pbmpi-master/data/pb_mpi -d nematode2.phy -T nematode2.phy.treefile  -nmax 15000 -poisson -cat -f nema1
