#!/bin/bash
#
#SBATCH --job-name=task1
#SBATCH --nodes=1    
#SBATCH --ntasks-per-node=8  
#SBATCH --time=00:20:00
#SBATCH --job-name=main
#SBATCH --output=out8.csv
#SBATCH --error=error8.txt

mpiexec fast.out 50 10000000000 0.1