#!/bin/bash
#
#SBATCH --nodes=1             # Number of nodes
#SBATCH --ntasks-per-node=8   # Number of processes per node
#SBATCH --job-name=MPI
# testing variables
 
echo "Запускаем тесты..."

# Очистка файла и добавление заголовка
echo "Processes N Time(s)" > out/logsSlug/resSched.txt

# Цикл по разным значениям N
my_it=1
N_array=(10 50 100)
for N in "${N_array[@]}"; do
    echo "Запуск с N = $N..." >> out/logsSlug/resSched.txt
    echo "===================================" >> out/logsSlug/resSched.txt
    
    # Цикл по количеству процессов
    for p in $(seq 1 8); do
        echo "  Запуск с $p процессами..."
        
	mpiexec -n $p ./gk.out $((2*N)) $N "/out/csv/${my_it}/{p}/" 0.004
            
    done
    my_it=$(($my_it + 1))

    echo "" >> out/logsSlug/resSched.txt  # Добавляем пустую строку для читабельности
done

echo "Тест завершен. Результаты записаны в resSched.txt" >> out/logsSlug/resSched.txt