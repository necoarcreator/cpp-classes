#!/bin/bash
#
#SBATCH --job-name=task1
#SBATCH --nodes=1    
#SBATCH --ntasks-per-node=8  
#SBATCH --time=00:20:00
#SBATCH --job-name=main
#SBATCH --output=out/logs/log8.txt
#SBATCH --error=out/logs/error8.txt
N=50
mpiexec gk.out $((2*N)) $N ./out/csv/test/8/ 0.01