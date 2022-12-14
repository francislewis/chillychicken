#!/bin/bash
# ================
# mpi_run.sh
# ================

#SBATCH --job-name=mpi_run
#SBATCH --partition=teach_cpu
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=28
#SBATCH --cpus-per-task=1
#SBATCH --time=0:30:0
#SBATCH --mem-per-cpu=250M
#SBATCH --array=1-56

# Load modules required for runtime e.g.
module load languages/intel/2020-u4

srun --mpi=pmi2 -n $SLURM_ARRAY_TASK_ID ./fft_mpi.o 32768 1
srun --mpi=pmi2 -n $SLURM_ARRAY_TASK_ID ./fft_mpi.o 32768 1
srun --mpi=pmi2 -n $SLURM_ARRAY_TASK_ID ./fft_mpi.o 32768 1
