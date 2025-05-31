#!/bin/bash
#
#SBATCH --job-name=task1
#SBATCH --nodes=1    
#SBATCH --ntasks-per-node=2  
#SBATCH --time=00:20:00
#SBATCH --job-name=main
#SBATCH --output=out2.txt
#SBATCH --error=error2.txt

mpiexec c.out 10 