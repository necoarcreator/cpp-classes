#include <stdio.h>
#include <iostream>
#include <omp.h>
#include <vector>
using namespace std;

int main(int argc, char** argv) {

	#ifdef _OPENMP
	printf("OpenMP is working correctly! %d\n", _OPENMP);
	#endif

	int i = 0;

	int myid, numprocs, numthreads;
	vector<int> vec = vector<int>(10); //get c-style array out of c++ vector, useful for mpi
	//cout << *vec.data() << endl;

	numprocs = omp_get_num_procs();
	cout << "Num of processors = " << numprocs << endl;

	omp_set_num_threads(10);

	numthreads = omp_get_num_threads();

	cout << "Num of threads = " << numthreads << endl;
	
	for (ptrdiff_t i = 0; i < 10; ++i) {
		vec[i] = i;
	}

	myid = omp_get_thread_num();
	
	cout << "consecutive part 1, myid = " << myid << endl;

	#pragma omp parallel shared(vec) private(myid, i) 
		{
			//директива для выполнения только мастером внутри параллельной части
			#pragma omp master 
				{
					printf("First val of vec = %d\n", vec[0]);
				}

			myid =  omp_get_thread_num();

			printf("Parallel part, myid = %d\n", myid);

		}

	myid = omp_get_thread_num();
	
	cout << "consecutive part 2, myid = " << myid << endl;


	return 0;

}