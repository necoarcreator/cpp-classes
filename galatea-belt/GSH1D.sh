#!/bin/bash
#
#SBATCH --job-name=GSH2D
#SBATCH --nodes=1
#SBATCH --ntasks=1             
#SBATCH --cpus-per-task=1 
#SBATCH --time=00:50:00

idx=$SLURM_ARRAY_TASK_ID


mpirun python3 gradsh_one.py
