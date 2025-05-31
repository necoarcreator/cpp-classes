#!/bin/bash
#
#SBATCH --nodes=1             # Number of nodes
#SBATCH --ntasks-per-node=8   # Number of processes per node
#SBATCH --job-name=MPI
# testing variables
 
echo "��������� �����..."

# ������� ����� � ���������� ���������
echo "Processes N Time(s)" > resSched.txt

# ���� �� ������ ��������� N
my_it=1
N_array=(2000 10000 50000)
for N in "${N_array[@]}"; do
    echo "������ � N = $N..." >> resSched.txt
    echo "===================================" >> resSched.txt
    
    # ���� �� ���������� ���������
    for p in $(seq 1 8); do
        echo "  ������ � $p ����������..."
        if [ $my_it -gt 1 ]
	then
        	
		mpiexec -n $p ./slow.out $N 1000000000 1e-4 >> "./res/${my_it}out_${p}.csv"
        else
        	mpiexec -n $p ./slow.out $N 1000000000 1e-4 >> "./res/out_${p}.csv"
	fi
    
    done
    my_it=$(($my_it + 1))

    echo "" >> resSched.txt  # ��������� ������ ������ ��� �������������
done

echo "���� ��������. ���������� �������� � resSched.txt" >> resSched.txt