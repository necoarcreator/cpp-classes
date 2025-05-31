#!/bin/bash
#
#SBATCH --nodes=1             # Number of nodes
#SBATCH --ntasks-per-node=8   # Number of processes per node
#SBATCH --job-name=MPI
# testing variables
 
echo "Запускаем тесты..."

# Очистка файла и добавление заголовка
echo "Processes N Time(s)" > resSched.txt

# Цикл по разным значениям N
my_it=1
N_array=(1000 1000000 100000000)
for N in "${N_array[@]}"; do
    echo "Запуск с N = $N..." >> resSched.txt
    echo "===================================" >> resSched.txt
    
    # Цикл по количеству процессов
    for p in $(seq 1 8); do
        echo "  Запуск с $p процессами..."
        if [ $my_it -gt 1 ]
	then
        	
		mpiexec -n $p ./c.out $N >> "./res/${my_it}out_${p}.txt"
        else
        	mpiexec -n $p ./c.out $N >> "./res/out_${p}.txt"
	fi
    
    done
    my_it=$(($my_it + 1))

    echo "" >> resSched.txt  # Добавляем пустую строку для читабельности
done

echo "Тест завершен. Результаты записаны в resSched.txt" >> resSched.txt