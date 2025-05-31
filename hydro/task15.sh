#!/bin/bash
#
#SBATCH --job-name=task1
#SBATCH --nodes=1    
#SBATCH --ntasks-per-node=5  
#SBATCH --time=00:20:00
#SBATCH --job-name=main
#SBATCH --output=out/logs/log5.txt
#SBATCH --error=out/logs/error5.txt
N=10
mpiexec weno.out $((4*N)) $N ./out/csv/test/5/"