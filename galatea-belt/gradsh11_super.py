#initialisation block
import numpy as np
import matplotlib.pyplot as plt
import time 
from mpi4py import MPI
import copy
import sys #for exiting


x_1 = -2
x_2 = 2
y_1 = -2
y_2 = 2
h_x, h_y, h = 0.01, 0.01, 0.01
tau = 0.00002
T = 4.2


N = int(T/tau)
x_grid = np.arange(x_1, x_2, h_x)
y_grid = np.arange(y_1, y_2, h_y)
I = len(x_grid)  
J = len(y_grid) 
X, Y = np.meshgrid(x_grid, y_grid)
PSI = np.zeros((J, I))
PSI_CH = np.zeros((I,J))
P = np.ones((I, J))
j_ex = np.zeros((I, J))
j_con = np.zeros((I, J))
H_x = np.zeros((I, J))
H_y = np.zeros((I, J))
j_all = np.zeros((I, J))
CHECK = None
end = None


p_0 = 1
r_c = 0.2
q = 0.2
j_0 = 2/r_c**2
σ = h**2

PSI_0 = 0

start_all = np.around(time.time())

for j in range(J): 
    for i in range (I): #считаем поэлементно ток проводников
        j_ex[j][i] = j_0 * (np.exp(-((x_grid[i] + 1)**2 + (y_grid[j])**2)/(r_c**2))         \
                        + np.exp(-((x_grid[i] - 1)**2 + (y_grid[j])**2)/(r_c**2)))
            
        PSI[j][i] = 4.8

if __name__ == "__main__":
    comm = MPI.COMM_WORLD # начинаем работу стандарта MPI
    nprocs = comm.Get_size() # индексы для каждого из потоков(процессоров)
    rank = comm.Get_rank() # число потоков
    ihproc = nprocs//2 + 1  # index_of_half_proc - так я коряво обозначил тот проц на который приходится элемент сетки I//2 и J//2
    sendbuf = None
    
    if rank == 0: # действия на корневом процессоре
       
        
        
        PSI = PSI.reshape((I * J, 1)) # раcпределяемый массив
        ave, res = divmod(PSI.size, nprocs) # определяем кол-во передаваемых элементов для каждого процессора
        
        count = [ave + 1 if p < res else ave for p in range(nprocs)] # составляем вектор из кол-ва отправляемых каждому процессору элементов
        
        count = np.array(count,dtype = int)
        #print(f'p = {rank} reports: \n\t{count}\nelements will be received by each process')

        displ = [sum(count[:p]) for p in range(nprocs)] # вектор индексов, начиная с которых будут передаваться данные
        displ = np.array(displ, dtype = int)
        PSI_std = - 1000;
        
    else:
        PSI = None
        count = np.zeros(nprocs, dtype = int) 
        displ = None
        
        #countCH = np.zeros(nprocs, dtype = int) 
        #displCH = None
        
    comm.Bcast(count, root = 0) # рассылаем из корневого процесса вектор кол-ва элементов
    
    
    
    CHECK =  np.zeros(nprocs)
    
    is_stab = 0

    recvbuf = np.zeros(count[rank]) # готовим пустой массив для приема данных; длина определяется соответствующим элементом массива count

    comm.Scatterv([PSI, count, displ, MPI.DOUBLE], recvbuf, root = 0) # рассылаем данные
    

    count = count//I # вектор из кол-ва отправляемых каждому процессору элементов теперь нумерует только строки
    
    for i in range(nprocs): #двумерные массивы нужной размерности для каждого процессора
        if rank == i:
            PSI = np.zeros((count[i], I)) # создаём двумерный массив обратно
            PSI_next = np.zeros((count[i], I)) 
            
   
    cumcount = np.cumsum(count)
    
    
    start = np.around(time.time())
    for n in range(N):
        if rank == ihproc: # для "срединного" проца вычисляем PSI_0 и отправляем на другие все процы
            PSI_0 = PSI[J//2 - cumcount[ihproc - 1],I//2]
            for k in range(nprocs):
                if k != ihproc: # на самого себя, понятное дело, не отправляем
                    comm.send(PSI_0, dest = k)
        else: # все остальные - получают
            PSI_0 = comm.recv(source = ihproc)
                    
       
            
        PSIleft = np.zeros(I) # левая мнимая строка 
        PSIright = np.zeros(I) # правая мнимая строка
           
        
        
        if rank == 0: #управляющий проц делится последней строкой и получает правую мнимую - необходмо для нашей разностной схемы
        
        
            for j in range(1, count[rank] - 1): #стандартные вычисления для всех строк всех процессоров
                for i in range (1, I - 1):
                    PSI_next[j][i] = PSI[j][i]             \
                                    + tau*((PSI[j][i + 1] - 2*PSI[j][i] + PSI[j][i - 1])/(h_x**2)             \
                                           + (PSI[j + 1][i] - 2*PSI[j][i] + PSI[j - 1][i])/(h_y**2)             \
                                               - 2* p_0/q**2 * (PSI[j][i] - PSI_0)*np.exp(-((PSI[j][i] - PSI_0)/q)**2)             \
                                                   + j_ex[j][i])
            comm.Send([PSI[count[rank] - 1], MPI.DOUBLE], dest = rank + 1) 
            comm.Recv([PSIright, MPI.DOUBLE], source = rank + 1) #помним, что для каждого проца за этой переменной кроются свои, разные, данные
            for i in range(1, I - 1): #отдельные вычисления для последней строки
                PSI_next[count[rank] - 1][i] = PSI[count[rank] - 1][i]             \
                                    + tau*((PSI[count[rank] - 1][i + 1] - 2*PSI[count[rank] - 1][i] + PSI[count[rank] - 1][i - 1])/(h_x**2)             \
                                           + (PSIright[i] - 2*PSI[count[rank] - 1][i] + PSI[count[rank] - 2][i])/(h_y**2)             \
                                               - 2* p_0/q**2 * (PSI[count[rank] - 1][i] - PSI_0)*np.exp(-((PSI[count[rank] - 1][i] - PSI_0)/q)**2)             \
                                                   + j_ex[count[rank] - 1][i])
            

        elif ((rank != 0)&(rank != nprocs - 1)): #для не первого и не последнего процев нужны мнимые строки с обеих сторон
            comm.Send([PSI[0], MPI.DOUBLE], dest = rank - 1) #отправляем первую строку на прош проц
            comm.Recv([PSIleft, MPI.DOUBLE], source = rank - 1) #получаем левую мнимую строку
            
            for i in range(1, I - 1): #вычисляем данные для первой строки
                PSI_next[0][i] = PSI[0][i]             \
                    + tau*((PSI[0][i + 1] - 2*PSI[0][i] + PSI[0][i - 1])/(h_x**2)             \
                           + (PSI[1][i] - 2*PSI[0][i] + PSIleft[i])/(h_y**2)             \
                               - 2* p_0/q**2 * (PSI[0][i] - PSI_0)*np.exp(-((PSI[0][i] - PSI_0)/q)**2)            \
                                   + j_ex[cumcount[rank - 1]][i])
            for j in range(1, count[rank] - 1): #стандартные вычисления для всех строк всех процессоров
                for i in range (1, I - 1):
                    PSI_next[j][i] = PSI[j][i]             \
                                    + tau*((PSI[j][i + 1] - 2*PSI[j][i] + PSI[j][i - 1])/(h_x**2)             \
                                           + (PSI[j + 1][i] - 2*PSI[j][i] + PSI[j - 1][i])/(h_y**2)             \
                                               - 2* p_0/q**2 * (PSI[j][i] - PSI_0)*np.exp(-((PSI[j][i] - PSI_0)/q)**2)             \
                                                   + j_ex[cumcount[rank - 1] + j][i])
            comm.Recv([PSIright, MPI.DOUBLE], source = rank + 1) #получаем правую мнимую строку
            comm.Send([PSI[count[rank] - 1], MPI.DOUBLE], dest = rank + 1) #отправляем последнюю строку на след проц
            for i in range(1, I - 1): # и для последней строки
                PSI_next[count[rank] - 1][i] = PSI[count[rank] - 1][i]             \
                    + tau*((PSI[count[rank] - 1][i + 1] - 2*PSI[count[rank] - 1][i] + PSI[count[rank] - 1][i - 1])/(h_x**2)             
                           + (PSIright[i] - 2*PSI[count[rank] - 1][i] + PSI[count[rank] - 2][i])/(h_y**2)             \
                               - 2* p_0/q**2 * (PSI[count[rank] - 1][i] - PSI_0)*np.exp(-((PSI[count[rank] - 1][i] - PSI_0)/q)**2)             \
                                   + j_ex[cumcount[rank - 1] + count[rank] - 1][i]) 
        else: #для последнего проца: делимся только первой строкой, получаем только левую мнимую
            comm.Send([PSI[0], MPI.DOUBLE], dest = rank - 1)
            comm.Recv([PSIleft, MPI.DOUBLE], source = rank - 1)
            for i in range(1, I - 1): #отдельно вычисляем данные для первой строки
                PSI_next[0][i] = PSI[0][i]             \
                    + tau*((PSI[0][i + 1] - 2*PSI[0][i] + PSI[0][i - 1])/(h_x**2)             \
                           + (PSI[1][i] - 2*PSI[0][i] + PSIleft[i])/(h_y**2)             \
                               - 2* p_0/q**2 * (PSI[0][i] - PSI_0)*np.exp(-((PSI[0][i] - PSI_0)/q)**2)             \
                                   + j_ex[cumcount[rank - 1]][i])
            for j in range(1, count[rank] - 1): #стандартные вычисления для всех строк всех процессоров
                for i in range (1, I - 1):
                    PSI_next[j][i] = PSI[j][i]             \
                                    + tau*((PSI[j][i + 1] - 2*PSI[j][i] + PSI[j][i - 1])/(h_x**2)             \
                                           + (PSI[j + 1][i] - 2*PSI[j][i] + PSI[j - 1][i])/(h_y**2)             \
                                               - 2* p_0/q**2 * (PSI[j][i] - PSI_0)*np.exp(-((PSI[j][i] - PSI_0)/q)**2)             \
                                               + j_ex[j + cumcount[rank - 1] - 1][i])
                 
            
                    
        PSI = copy.deepcopy(PSI_next) # переход на след шаг по времени, если установления нет
        
        
              
    recvbuf2 = np.zeros((J, I)) # готовимся к приему на корневом процессоре результатов
    #print(f'recvbuf2 is {recvbuf2}, the size is {len(recvbuf2)}')
    
    count = count * I # вектор из кол-ва отправляемых каждому процессору элементов  
    PSI_next = PSI_next.reshape((count[rank], 1))
    comm.Gatherv([PSI_next, count[rank], MPI.DOUBLE], [recvbuf2, count, displ, MPI.DOUBLE], root=0) #даем всем команду на рассылку и сбор данных корневым процессом
    PSI_new = recvbuf2.reshape((J, I))
    end = np. around(time.time() - start, 3)/60
    if rank == 0:

        
        with open('/home/saa128/galatea-belt/properties/exectime.txt', 'a') as f:
            f.write(f"Maximum pressure is p_0 = {p_0} units. The calculation time was about t = {end} minutes with T = {T} units, the number of threads was n = {nprocs}.")  #сохраняем время вычисления
            np.savetxt(r'/home/saa128/galatea-belt/properties/PSI_new.txt', PSI_next)  #сохраняем на всякий пожарный получившийся массив пси
        
      
    if rank == 0:

        
        #visualisation block


        col = ['p_0', 'h', 'tau', 'T', 'σ', 'r_c', 'q']
        width = [0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2]
        row = ['values']
        cells = [[str(p_0), str(h), str(tau), str(T), str(np.around(σ,4)), str(r_c), str(q)]]

        manual_locations = [(-1.5,0.), (-0.75,0.1), (-1,-0.2), 
                            (1.5, 0.), (0.75, 0.1), (1.,-0.2)]
        CS = plt.contour(X, Y, PSI_new)
        plt.clabel(CS, manual = manual_locations)
        plt.title("Grad-Shafranov equation for axial symmetric galatea", fontsize = 16)
        the_table = plt.table(rowLabels=row,
                              colLabels=col,
                              cellText = cells,
                              colWidths = width,
                              colLoc = 'center',
                              cellLoc = 'center',
                              loc = 'upper left',
                              bbox = (1.2 , 0.3, 0.9, 0.6),) 
        plt.savefig('/home/saa128/galatea-belt/properties/G-Sh_field')
        plt.show()


        for i in range (1, I - 1):
            for j in range (1, J - 1):
                P[j][i] = p_0 * np.exp(-((PSI_new[j][i] - PSI_0)**2)/q**2)
                H_x[j][i] = (PSI_new[j + 1][i] - PSI_new[j - 1][i])/(2*h_x)
                H_y[j][i] = - (PSI_new[j][i + 1] - PSI_new[j][i - 1])/(2*h_y)
                
        for i in range (1, I - 1):
            H_x[0][i] = (PSI_new[1][i] - PSI_new[0][i])/h_x
            H_y[0][i] = - (PSI_new[0][i + 1] - PSI_new[0][i])/h_y
        for j in range (1, J - 1):
            H_x[j][0] = (PSI_new[j + 1][0] - PSI_new[j][0])/h_x
            H_y[j][0] = - (PSI_new[j][1] - PSI_new[j][0])/h_y
        for j in range (1, J - 1):        
            H_x[j][I - 1] = (PSI_new[j][I - 1] - PSI_new[j - 1][I - 1])/h_x
            H_y[j][I - 1] = - (PSI_new[j][I - 1] - PSI_new[j][I - 2])/h_y
        for i in range (1, I - 1):        
            H_x[J - 1][i] = (PSI_new[J - 1][i] - PSI_new[J - 2][i])/h_x
            H_y[J - 1][i] = - (PSI_new[J - 1][i] - PSI_new[J - 1][i - 1])/h_y

        H_x *= -1
        H_y *= -1
        manual_locations = [(0.,1.6), (0,-1.8), (-1,-0.7), 
                            (0., 0.), (1.,-0.2), (1.55,0.2)]
        fig, ax = plt.subplots(1, 2,figsize=(15,5))
        KC = ax[0].contour(X, Y, P)
        ax[0].clabel(KC, manual = manual_locations)
        ax[0].set_title('A visualistion for pressure levels', fontsize = 16)
        ax[1].imshow(P)
        the_table = plt.table(rowLabels=row,
                              colLabels=col,
                              cellText = cells,
                              colWidths = width,
                              colLoc = 'center',
                              cellLoc = 'center',
                              loc = 'upper left',
                              bbox = (1.4 , 0.8, 0.6, 0.2),) 
        ax[1].set_title('A visualisation for pressure split', fontsize = 16)
        plt.savefig('/home/saa128/galatea-belt/properties/Pressure')
        p2 = ax[1].imshow(P, cmap='plasma', aspect='equal', interpolation='gaussian', origin="lower", extent=(0, 80, 0, 80))
        fig.colorbar(p2, ax=ax[1])


        fig, ax = plt.subplots()
        ax.quiver(H_x, H_y, color = 'red')
        fig.set_figwidth(10)   
        fig.set_figheight(10)
        ax.set_title('A configuration of magnetic field', fontsize = 16)
        the_table = plt.table(rowLabels=row,
                              colLabels=col,
                              cellText = cells,
                              colWidths = width,
                              colLoc = 'center',
                              cellLoc = 'center',
                              loc = 'upper left',
                             bbox = (1.1 , 0.8, 0.6, 0.2),) 
        plt.savefig('/home/saa128/galatea-belt/properties/Vctr_mgntc_fld')

        plt.show()


        fig, ax = plt.subplots(1, 2,figsize=(15,5))
        manual_locations = [(1., 0.), (0,-1.5), (-1,-0.6), 
                            (-1., 0.3), (1.,-0.4), (1.6,0.)]
        KC = ax[0].contour(X,Y, PSI_new - PSI_0)
        ax[0].clabel(KC, manual = manual_locations)
        BC = ax[1].contour(X, Y, P)
        manual_locations = [(0., 0.5), (0,-1.5), (-1,-0.6), 
                            (-1., 0.3), (0.,-0.4), (1.6,0.)]
        ax[1].clabel(BC, manual = manual_locations)
        ax[0].set_title("Levels of magnetic flux", fontsize = 16)
        ax[1].set_title("Plasma pressure", fontsize = 16)
        ax[0].legend(loc=3,bbox_to_anchor=(0.1, 1.1, 0.5, 0.5), fontsize = 16)
        the_table = plt.table(rowLabels=row,
                              colLabels=col,
                              cellText = cells,
                              colWidths = width,
                              colLoc = 'center',
                              cellLoc = 'center',
                              loc = 'upper left',
                              bbox = (1.1 , 0.8, 0.6, 0.2),) 
        plt.savefig('/home/saa128/galatea-belt/properties/pic_1')
        plt.show()


        for i in range (1, I - 1):
            for j in range (1, J - 1):
                j_all[j][i] = (H_y[j][i + 1] - H_y[j][i - 1])/(2*h)                   \
                            - (H_x[j + 1][i] - H_x[j - 1][i])/(2*h) 
        j_pl = j_all - j_ex




        fig, ax = plt.subplots(1, 2,figsize=(15,5))
        manual_locations = [(0., 0.), (0, -1.7), (-1, -0.4), 
                            (-1., 0.3), (1., -0.2), (1., 0.1)]
        ax[0].quiver(H_x, H_y, color = 'red')
        DC = ax[1].contour(X, Y, j_pl)
        ax[1].clabel(DC, manual = manual_locations)
        ax[0].set_title("H-vector in plasma", fontsize = 16)
        ax[1].set_title("Plasma current density", fontsize = 16)
        the_table = plt.table(rowLabels=row,
                              colLabels=col,
                              cellText = cells,
                              colWidths = width,
                              colLoc = 'center',
                              cellLoc = 'center',
                              loc = 'upper left',
                              bbox = (1.1 , 0.8, 0.6, 0.2),) 
        plt.savefig('/home/saa128/galatea-belt/properties/pic_2')
    
    end_all = np. around(time.time() - start_all, 3)/60

    if (rank == 0):
        with open('/home/saa128/galatea-belt/properties/execsr.txt', 'a') as f:
            f.write(f"The parallel calculation time was about t = {end} minutes, the overall calculation time was {end_all} minutes. The proportion of parallel computations was p = {end/end_all}. The number of threads was n = {nprocs}.")  #сохраняем время вычисления
 

