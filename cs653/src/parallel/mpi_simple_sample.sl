#!/bin/bash
#SBATCH --ntasks-per-node=2
#SBATCH --nodes=1
#SBATCH --time=00:00:10
#SBATCH --output=mpi_simple.out
#SBATCH -A lc_an2
WORK_HOME=/home/rcf-proj/an2/yourID
cd $WORK_HOME
srun -n $SLURM_NTASKS --mpi=pmi2 ./mpi_simple
