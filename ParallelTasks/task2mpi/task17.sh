#!/bin/bash
#
#SBATCH --job-name=task1
#SBATCH --nodes=1    
#SBATCH --ntasks-per-node=7  
#SBATCH --time=00:20:00
#SBATCH --job-name=main
#SBATCH --output=out7.csv
#SBATCH --error=error7.txt

mpiexec fast.out 2000 1000000 0.0001