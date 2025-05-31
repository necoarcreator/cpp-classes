#!/bin/bash
#
#SBATCH --nodes=1             # Number of nodes
#SBATCH --ntasks-per-node=8   # Number of processes per node
#SBATCH --job-name=MPI
# testing variables
 
echo "��������� �����..."

# ������� ����� � ���������� ���������
echo "Processes N Time(s)" > out/logsSlug/resSched.txt

# ���� �� ������ ��������� N
my_it=1
N_array=(10 50 100)
for N in "${N_array[@]}"; do
    echo "������ � N = $N..." >> out/logsSlug/resSched.txt
    echo "===================================" >> out/logsSlug/resSched.txt
    
    # ���� �� ���������� ���������
    for p in $(seq 1 8); do
        echo "  ������ � $p ����������..."
        
	mpiexec -n $p ./gk.out $((2*N)) $N "/out/csv/${my_it}/{p}/" 0.004
            
    done
    my_it=$(($my_it + 1))

    echo "" >> out/logsSlug/resSched.txt  # ��������� ������ ������ ��� �������������
done

echo "���� ��������. ���������� �������� � resSched.txt" >> out/logsSlug/resSched.txt