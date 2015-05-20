#include <stdio.h>
#include <math.h>
#include <mpi.h>

int main(int argc, char *argv[])
{
	MPI_Init(NULL, NULL);

	int array_size = 1000000;
	int iterations = 1000;

	int rank, size;

	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &size);

	MPI_Finalize();

	return;
}
