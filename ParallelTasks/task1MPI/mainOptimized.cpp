#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char** argv) {
    int size, rank;
    
    long N = atol(argv[1]);
    double L = 1.0;
    double summ = 0.0;
    double dx = L / N;  

    double borderL, borderR;  
    MPI_Status Status;
    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // начальное время
    MPI_Barrier(MPI_COMM_WORLD);
    double start = MPI_Wtime();

    // размер участка 
    int chunk = N / size + (rank < N % size ? 1 : 0);

    if (rank == 0) {
        // 0 рассчитывает границы 
        for (int p = 0; p < size; ++p) {
            double left = (p * (N / size) + (p < N % size ? p : N % size)) * dx;
            double right = left + (chunk) * dx;

            if (p == size - 1) right = L;  

            if (p == 0) {
                borderL = left;
                borderR = right;
            } else {
                MPI_Send(&left, 1, MPI_DOUBLE, p, 0, MPI_COMM_WORLD);
                MPI_Send(&right, 1, MPI_DOUBLE, p, 1, MPI_COMM_WORLD);
            }
        }
    } else {
       
        MPI_Recv(&borderL, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &Status);
        MPI_Recv(&borderR, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &Status);
    }

    // вычисляем интеграл 
    double x = borderL;
    for (int k = 0; k < chunk; ++k) {
        double f1 = 4.0 / (1 + x * x);
        double f2 = 4.0 / (1 + (x + dx) * (x + dx));
        summ += 0.5 * (f1 + f2) * dx;
        x += dx;
    }

    
    double total_sum = 0.0;
    MPI_Reduce(&summ, &total_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    
    if (rank == 0) {
        printf("The summ is %f\n", total_sum);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double end = MPI_Wtime();
    printf("%d, %f\n", rank, end - start);

    MPI_Finalize();
    return 0;
}