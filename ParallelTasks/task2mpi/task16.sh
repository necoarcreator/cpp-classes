#!/bin/bash
#
#SBATCH --job-name=task1
#SBATCH --nodes=1    
#SBATCH --ntasks-per-node=6  
#SBATCH --time=00:20:00
#SBATCH --job-name=main
#SBATCH --output=out6.csv
#SBATCH --error=error6.txt

mpiexec slow.out 50 1000000 0.1