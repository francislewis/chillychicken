#!/bin/bash
# ================
# omp_run.sh
# ================

#SBATCH --job-name=omp_run
#SBATCH --partition=teach_cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28
#SBATCH --time=0:30:0
#SBATCH --mem=300M
#SBATCH --array=1-28

# Load modules required for runtime e.g.
module load languages/intel/2020-u4

srun ./fft_omp.o 32768 1 $SLURM_ARRAY_TASK_ID
srun ./fft_omp.o 32768 1 $SLURM_ARRAY_TASK_ID
srun ./fft_omp.o 32768 1 $SLURM_ARRAY_TASK_ID
