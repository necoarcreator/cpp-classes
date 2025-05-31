#!/bin/bash
#
#SBATCH --job-name=task1
#SBATCH --nodes=1    
#SBATCH --ntasks-per-node=1  
#SBATCH --time=00:20:00
#SBATCH --job-name=main
#SBATCH --output=out/logs/log1.txt
#SBATCH --error=out/logs/error1.txt

N=50
mpiexec gk.out $((2*N)) $N ./out/csv/test/1/ 0.05

