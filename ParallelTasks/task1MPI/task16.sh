#!/bin/bash
#
#SBATCH --job-name=task1
#SBATCH --nodes=1    
#SBATCH --ntasks-per-node=6  
#SBATCH --time=00:20:00
#SBATCH --job-name=main
#SBATCH --output=out6.txt
#SBATCH --error=error6.txt

mpiexec c.out 100000000