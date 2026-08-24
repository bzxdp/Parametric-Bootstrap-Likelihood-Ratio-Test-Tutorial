#!/bin/bash
#SBATCH --job-name=CAT_pmsf
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=60GB
#SBATCH --time=02:00:00


module load gnu12/12.2.0
module load openmpi4/4.1.5

export OMPI_MCA_pml=ob1
export OMPI_MCA_btl=self,tcp


cd $SLURM_SUBMIT_DIR


srun -n 16 /projects/b56o/software/PATCHED_PB_MPI_MICAT/pbmpi/data/readpb_mpi -x 6000 67 12755 -r nema2


