#!/bin/bash
# ================
# serial_run.sh
# ================

#SBATCH --job-name=serial_run
#SBATCH --partition=teach_cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:30:0
#SBATCH --mem=500M
#SBATCH --array=2-131072

# Load modules required for runtime e.g.
module load languages/intel/2020-u4

srun ./fft_serial.o $SLURM_ARRAY_TASK_ID 2
srun ./fft_serial.o $SLURM_ARRAY_TASK_ID 2
srun ./fft_serial.o $SLURM_ARRAY_TASK_ID 2
