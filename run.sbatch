#!/bin/bash
#SBATCH --job-name=4_6E-5_9_0.2_700
#SBATCH --output=MITgcm_sbatch.out
#SBATCH --error=MITgcm_sbatch.err
#SBATCH --time=32:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --exclusive

module unload intel intelmpi
module unload netcdf
module load netcdf/4.4.1+intel-16.0
module load netcdf-fortran/4.4.4+intel-16.0
module load intelmpi/5.1+intel-16.0

mpirun -n 28 ./mitgcmuv

