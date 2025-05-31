#!/bin/bash
#
#SBATCH --job-name=task1
#SBATCH --nodes=1    
#SBATCH --ntasks-per-node=4  
#SBATCH --time=00:20:00
#SBATCH --job-name=main
#SBATCH --output=out4.csv
#SBATCH --error=error4.txt

mpiexec fast.out 100 1000000 0.1