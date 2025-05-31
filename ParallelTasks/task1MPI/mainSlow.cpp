#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
 
//#include <vector>

using namespace std;

int main(int argc, char** argv) {

	int size, rank;
	//берЄм число точек из аргументов командной строки
	long N = atol(argv[1]);
	double L = 1.0;
	double summ = 0.0;
	double dx = L / N;

	double* xp;
	double* x_locp;
	double x_right;


	MPI_Status Status;
	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	//замер€ем начальное врем€
	MPI_Barrier(MPI_COMM_WORLD);
	double start = MPI_Wtime();
	//инициализируем массив x и локальные массивы x_loc, запол€нем x в управл€ющем процессе, рассылаем его в x_loc
	//рассылаем свои части x_loc из x каждому процессу, считаем интеграл, отправл€ем результат на управл€ющий процесс 
	

	int* sendcounts = (int*) calloc(size, sizeof(int));
    	int* displs = (int*) calloc(size, sizeof(int));

	
	size_t offset = 0;
	for (int i = 0; i < size; i++) {
        	sendcounts[i] = N / size + (i < N % size ? 1 : 0);
        	displs[i] = offset;
        	offset += sendcounts[i];
    	}

	int chunk = sendcounts[rank];
	x_locp = (double*) calloc(chunk, sizeof(double));

	if (rank == 0) {

		xp = (double*)calloc(N, sizeof(double));


		for (ptrdiff_t j = 0; j < N; ++j) {
			xp[j] = j * dx;
		}

		for (ptrdiff_t p = 1; p < size; ++p) {
			MPI_Send(&xp[displs[p]], sendcounts[p], MPI_DOUBLE, p, p, MPI_COMM_WORLD);
		}
		for (ptrdiff_t j = 0; j < chunk; ++j) {
			x_locp [j] = xp[j];
		}
	}
	else {
		MPI_Recv(x_locp, chunk, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD, &Status);
	}

	
	if (rank % 2 == 0) {  //  чЄтные сначала принимают 
		
		if (rank != size - 1) {
		//printf("recieving from %d to %d\n", rank + 1, rank);
        	MPI_Recv(&x_right, 1, MPI_DOUBLE, rank + 1, MPI_ANY_TAG, MPI_COMM_WORLD, &Status);
		//printf("recieved from %d to %d\n", rank + 1, rank);
    		}
		else {
    		x_right = N * dx;
		}
	} 
	else if (rank != 0) {  // нечЄтные сначала отправл€ют
		
		//printf("sending from %d to %d\n", rank, rank - 1); 
		MPI_Send(&x_locp[0], 1, MPI_DOUBLE, rank - 1, rank - 1, MPI_COMM_WORLD);
		//printf("sended from %d to %d\n", rank, rank - 1);
	}
	
	//printf("=======NEW TACT=======\n");
	
	if (rank % 2 != 0) {  // нечЄтные принимают 
		
		if (rank != size - 1) {
        	//printf("recieving from %d to %d\n", rank + 1, rank);
        	MPI_Recv(&x_right, 1, MPI_DOUBLE, rank + 1, MPI_ANY_TAG, MPI_COMM_WORLD, &Status);
		//printf("recieved from %d to %d\n", rank + 1, rank);
    		}
		else {
    		x_right = N * dx;
		}
	} 
	else if (rank != 0) {  // чЄтные отправл€ют 	 	
	
		//printf("sending from %d to %d\n", rank, rank - 1); 
		MPI_Send(&x_locp[0], 1, MPI_DOUBLE, rank - 1, rank - 1, MPI_COMM_WORLD);
		//printf("sended from %d to %d\n", rank, rank - 1);
	}


	//printf("x_0 = %f and x_n = %f\n", x_locp[0], x_locp[chunk - 1]);

	for (ptrdiff_t k = 0; k < chunk - 1; ++k) {
		summ += 0.5 * (4.0 / (1 + x_locp[k] * x_locp[k]) + 4.0 / (1 + x_locp[k+1] * x_locp[k+1])) * dx;
	}
	summ += 0.5 * (4.0 / (1 + x_locp[chunk - 1] * x_locp[chunk - 1]) + 4.0 / (1 + x_right * x_right)) * dx;
	
	if (rank != 0) {
		MPI_Send(&summ, 1, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD);
	}
	
	//принимаем результаты, суммируем, выводим
	if (rank == 0) {

		double sum_temp = 0.0;

		for (size_t p = 1; p < size; ++p) {
			MPI_Recv(&sum_temp, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &Status);
			summ += sum_temp;
		}

		printf("The summ is %f", summ);
		printf("\n");
	}

	free(x_locp);
	free(sendcounts);
	free(displs);

	MPI_Barrier(MPI_COMM_WORLD);
	double end = MPI_Wtime();
	printf("%d, %f", rank, end - start);
	printf("\n");

	MPI_Finalize();
	return 0;
}