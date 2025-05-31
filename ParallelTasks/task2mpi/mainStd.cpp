#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <iostream>
#include <math.h>
#define _USE_MATH_DEFINES

#define _DEBUG_

using namespace std;

int main(int argc, char** argv) {
	//температуропроводность
	double k = 1.0;

	int size, rank;
	//берём число точек из аргументов командной строки
	long N = stol(argv[1]);
	double L = 1.0, T = 1.0, T_stop = stod(argv[3]);
	double u_0 = 1.0;
	bool isT10 = false;
	//шаг по времени по CFL = 1/2
	double dx = L / N, dt = T * dx * dx / (2 * k);
	double borderL, borderR;
	//фиктивные ячейки
	double u_left, u_right;

	MPI_Status Status;
	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	double* u_theor = new double[11];
	for (int i = 0; i < 10; ++i) u_theor[i] = 0.0;

	
	for (ptrdiff_t i = 0; i < 11; ++i) {

		double x = i * (L / 10);
		for (ptrdiff_t m = 0; m < 2; ++m) {
			
		u_theor[i] += 4 * u_0 / M_PI * exp(- k * M_PI * M_PI * (2 * m + 1) * (2 * m + 1) *
												T_stop / L / L) / (2 * m + 1) * sin(M_PI * (2 * m + 1) * x / L);
		}
		
		x += L / 10;
	}
	//замеряем начальное время
	MPI_Barrier(MPI_COMM_WORLD);
	double start = MPI_Wtime(); 

	//размер разделения
	int chunk = (N + 1) / size + (rank < (N + 1) % size ? 1 : 0);
	double* up = new double[chunk];
	double* u_newp = new double[chunk];
	
	if (rank == 0) {

		for (int p = 0; p < size; ++p) {
			double left = (p * ((N + 1) / size) + (p < (N + 1) % size ? p : (N + 1) % size)) * dx;

			double right = left + ((N + 1) / size + (p < (N + 1) % size ? 1 : 0))*dx;

			if (p == size - 1) right = L;

			if (p == 0) {
				borderL = left;
				borderR = right;
			}
			else {
				MPI_Send(&left, 1, MPI_DOUBLE, p, 0, MPI_COMM_WORLD);
				MPI_Send(&right, 1, MPI_DOUBLE, p, 1, MPI_COMM_WORLD);
			}
		}
	}
	else {

		MPI_Recv(&borderL, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &Status);
		MPI_Recv(&borderR, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &Status);

	}

	double time = 0.0;
	size_t iter = 0, max_iter = stoll(argv[2]);

	//начальные условия
	for (ptrdiff_t j = 0; j < chunk; ++j) {
		up[j] = u_0;
		u_left = u_0;
		u_right = u_0;
		u_newp[j] = 0.0;

	}
	

	while (time < 2 * T_stop && iter < max_iter) {

		//граничные условия и прямой проход

		//O(p)
		if (rank > 0) {
			MPI_Send(&up[0], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
		}
		if (rank < size - 1) {
			MPI_Recv(&u_right, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &Status);
		}
		if (rank < size - 1) {
			MPI_Send(&up[chunk - 1], 1, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD);
		}
		if (rank > 0) {
			MPI_Recv(&u_left, 1, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, &Status);
		}

		if (rank == 0) {
			u_left = 0.0;
			up[0] = 0.0;
			//if (iter % 10000 == 0) { 
			//	cout << time << ", " << iter << endl;
			//}
		}
		if (rank == size - 1) {
			u_right = 0.0;
			up[chunk - 1] = 0.0;
		}

		u_newp[0] = up[0] + k * dt / (dx * dx) * (up[1] - 2 * up[0] + u_left);

		for (ptrdiff_t j = 1; j < chunk - 1; ++j) {

			u_newp[j] = up[j] + k * dt / (dx * dx) * (up[j + 1] - 2 * up[j] + up[j - 1]);
		}

		u_newp[chunk - 1] = up[chunk - 1] + k * dt / (dx * dx) * (u_right - 2 * up[chunk - 1] + up[chunk - 2]);

		//обновление значений переменной
		for (ptrdiff_t j = 0; j < chunk; ++j) {

			up[j] = u_newp[j];
		}

		time += dt;
		++iter;
		if (time - T_stop > 0.0 and not isT10) {
			isT10 = true;
			MPI_Barrier(MPI_COMM_WORLD);
			double end = MPI_Wtime();
			double* u_finp = nullptr;
			int* recvcounts = nullptr;
			int* displs = nullptr;
			int offset = 0;
			
			if (rank == 0)  {
			
				u_finp = new double[N + 1];
				recvcounts = new int[size];
				displs = new int[size];

				for (int i = 0; i < size; ++i) {
					recvcounts[i] = (N + 1) / size + (i < (N + 1) % size ? 1 : 0);
					displs[i] = offset;
					offset += recvcounts[i];
				}

			}
		
			MPI_Gatherv(up, chunk, MPI_DOUBLE,
				u_finp, recvcounts, displs, MPI_DOUBLE,
				0, MPI_COMM_WORLD);

			if (rank == 0) {

					double x = 0.0;
					cout  << size << ", " << end - start << endl;
					cout  << "x, " << "u, " << "u_an" << endl;

					for (ptrdiff_t i = 0; i < 11; ++i) {

						cout << x << ", " << u_finp[i * N / 10] << ", " << u_theor[i] << endl;
						x += L / 10;

					}			
			delete[] u_finp;
			delete[] recvcounts;
			delete[] displs; 
			}

		#ifdef _DEBUG_
		iter = max_iter;
		#endif
		}
	}

		delete[] up;
		delete[] u_newp;
		delete[] u_theor;

	
	MPI_Finalize();
	return 0;
}