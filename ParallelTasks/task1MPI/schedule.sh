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
N_array=(1000 1000000 100000000)
for N in "${N_array[@]}"; do
    echo "������ � N = $N..." >> resSched.txt
    echo "===================================" >> resSched.txt
    
    # ���� �� ���������� ���������
    for p in $(seq 1 8); do
        echo "  ������ � $p ����������..."
        if [ $my_it -gt 1 ]
	then
        	
		mpiexec -n $p ./c.out $N >> "./res/${my_it}out_${p}.txt"
        else
        	mpiexec -n $p ./c.out $N >> "./res/out_${p}.txt"
	fi
    
    done
    my_it=$(($my_it + 1))

    echo "" >> resSched.txt  # ��������� ������ ������ ��� �������������
done

echo "���� ��������. ���������� �������� � resSched.txt" >> resSched.txt