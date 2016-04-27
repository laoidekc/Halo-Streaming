#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

void calculation_delay(double*, int);

//This code is for testing the injection and data retrieval rates of MPI
int main(int argc, char *argv[])
{
	MPI_Init(NULL,NULL);
	
	int rank, size;
	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &size);

	if(size != 2 && rank == 0)
	{
		printf("Must be run on exactly two processes\n");
	}

	if(argc < 3 && rank == 0)
	{
		printf("Need to supply number of messages, message size, and amount of calculation as arguments\n");
		return 0;
	}

	int num_messages = atoi(argv[0]);
	int message_size = atoi(argv[1]);
	int calculation_size = atoi(argv[2]);

	int i;

	double message_data[message_size];

	// Send tests

	double times[num_messages];
	double calculation_data[calculation_size];
	for(i=0;i<calculation_size;i++)
	{
		calculation_data[i] = rand()%10000;
	}

	if(rank==0)
	{
		times[0]=MPI_Wtime();
		for(i=0;i<num_messages;i++)
		{
			MPI_Send(message_data,message_size,MPI_DOUBLE,1,i,comm);
			calculation_delay(calculation_data,calculation_size);
			times[i+1]=MPI_Wtime();
		}
	}
	// Receive tests
}

void calculation_delay(double *calculation_data, int calculation_size)
{
	for(int i=0;i<calculation_size;i++)
	{
		calculation_data[i] = calculation_data[i] + 1;
	}
}
