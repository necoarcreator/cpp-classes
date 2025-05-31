#!/bin/bash
#
#SBATCH --job-name=task1
#SBATCH --nodes=1    
#SBATCH --ntasks-per-node=3  
#SBATCH --time=00:20:00
#SBATCH --job-name=main
#SBATCH --output=out3.txt
#SBATCH --error=error3.txt

mpiexec c.out 10 