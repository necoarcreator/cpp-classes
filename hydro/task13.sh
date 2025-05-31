#!/bin/bash
#
#SBATCH --job-name=task1
#SBATCH --nodes=1    
#SBATCH --ntasks-per-node=3  
#SBATCH --time=00:20:00
#SBATCH --job-name=main
#SBATCH --output=out/logs/log3.txt
#SBATCH --error=out/logs/error3.txt
N=10
mpiexec weno.out $((4*N)) $N ./out/csv/test/3/"