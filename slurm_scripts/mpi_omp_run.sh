#!/bin/bash
# ================
# mpi_omp_run.sh
# ================

#SBATCH --job-name=mpi_omp_run
#SBATCH --partition=teach_cpu
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=7
#SBATCH --time=0:30:0
#SBATCH --mem-per-cpu=250M

# Load modules required for runtime e.g.
module load languages/intel/2020-u4

srun --mpi=pmi2 -n 4 ./fft_mpi_omp.o 32768 1 7
srun --mpi=pmi2 -n 4 ./fft_mpi_omp.o 32768 1 7
srun --mpi=pmi2 -n 4 ./fft_mpi_omp.o 32768 1 7
