#!/bin/bash
#
#SBATCH --job-name=GSH2D
#SBATCH --nodes=2
#SBATCH --ntasks=25            
#SBATCH --cpus-per-task=1 
#SBATCH --time=00:30:00

idx=$SLURM_ARRAY_TASK_ID


mpirun python3 gradsh11_super.py
