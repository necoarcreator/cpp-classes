#!/bin/bash
#
#SBATCH --job-name=task1
#SBATCH --nodes=1    
#SBATCH --ntasks-per-node=1  
#SBATCH --time=00:20:00
#SBATCH --job-name=main
#SBATCH --output=out1.csv
#SBATCH --error=error1.txt

mpiexec fast.out 50 1000000 0.1