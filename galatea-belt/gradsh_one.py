#initialisation block
import numpy as np
import time 
from mpi4py import MPI

#geometry block
x_1 = -2
x_2 = 2
y_1 = -2
y_2 = 2
h_x = h_y = h = 0.05
tau = 0.0001
T = 4.2

#grid block
N = int(T/tau)
x_grid = np.arange(x_1, x_2 + h_x, h_x)
y_grid = np.arange(y_1, y_2 + h_y, h_y)
I = len(x_grid)  
J = len(y_grid) 
X, Y = np.meshgrid(x_grid, y_grid)
PSI = np.zeros((I, J))
PSI_next = np.zeros((I, J))
j_ex = np.zeros((I, J))

#parameters block
p_0 = 1.0
r_c = 0.2
q = 0.2
j_0 = 2/r_c**2
for i in range (1, I - 1):
        for j in range (1, J - 1):
            j_ex[j][i] = j_0* \
            (np.exp(-((x_grid[i] + 1)**2 + (y_grid[j])**2)/(r_c**2)) \
            + np.exp(-((x_grid[i] - 1)**2 + (y_grid[j])**2)/(r_c**2)))

#numeric experiment block

if __name__ == "__main__":
    comm = MPI.COMM_WORLD # начинаем работу стандарта MPI
    nprocs = comm.Get_size() # индексы для каждого из потоков(процессоров)
    rank = comm.Get_rank()
    
    start = time.time()
    for n in range(N + 1):
        PSI_0 = PSI[J//2][I//2] 
        if not n%1000: 
            print(f'The step is {int(n/1000)}, magnetic flux in centre equals to {np.around(PSI_0, 3)}, {np.around(time.time() - start)/60} minutes have already passed')
            print(np.max(PSI_next))
        for i in range (1, I - 1):
            for j in range (1, J - 1):
                PSI_next[j][i] = PSI[j][i] \
                    + tau*((PSI[j][i + 1] - 2*PSI[j][i] + PSI[j][i - 1])/(h_x**2) \
                    + (PSI[j + 1][i] - 2*PSI[j][i] + PSI[j - 1][i])/(h_y**2) \
                    - 2* p_0/q**2 * (PSI[j][i] - PSI_0)*np.exp(-((PSI[j][i] - PSI_0)/q)**2) \
                    + j_ex[j][i])
        if n != N:    
            PSI = PSI_next.copy()
    print('block was initialised')
    end = np.around((time.time() - start)/60,3) 
    with open('exectime.txt', 'w') as f:
        f.write(f"Maximum pressure is p_0 = {p_0} units. The calculation time was about t = {end} minutes with T = {T} units, the number of threads was n = 1.")  #сохраняем время вычисления
        np.savetxt(r'PSI_new.txt', PSI_next)  #сохраняем на всякий пожарный получившийся массив пси

