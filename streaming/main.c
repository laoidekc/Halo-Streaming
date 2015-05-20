#include <stdio.h>
#include <math.h>
#include <mpi.h>

int main(int argc, char *argv[])
{
	MPI_Init(NULL, NULL);

	int array_size = 1000000;
	int iterations = 1000;
	int send_buffer_size = 10;

	int rank, size, i, j, location = array_size-2, tag;
	double local_data[array_size+2], start_time, time, total = 0, checksum;

	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &size);
	MPI_Request request[send_buffer_size];
	MPI_Status status[send_buffer_size];

	int neighbour_left = (rank + size - 1) % size;
	int neighbour_right = (rank + 1) % size;

	srand(rank);

	for(i=0;i<array_size;i++)
	{
		local_data[i] = (double)i;
	}

	for(tag=0;tag<send_buffer_size;tag++)
	{
		MPI_Send_init(local_data,2,MPI_DOUBLE,neighbour_left,tag,comm,&request[tag]);	
	}

	tag = 0;

	if(rank == 0)
	{
		start_time = MPI_Wtime();
	}

	for(j=0;j<iterations;j++)
	{
		MPI_Wait(&request[tag],&status[tag]);
		MPI_Start(&request[tag]);
		tag = (tag + 1) % send_buffer_size;

		for(i=0;i<location;i++)
		{
			local_data[i] = (local_data[i] + local_data[i+1] + local_data[i+2])/3;
		}
		location -= 2;
	}

	location = array_size - 2;
	tag = 0;

	for(j=0;j<iterations;j++)
	{
		MPI_Recv(&local_data[array_size],2,MPI_DOUBLE,neighbour_right,tag,comm,&status[tag]);
		tag = (tag + 1) % send_buffer_size;

		for(i=location;i<array_size;i++)
		{
			local_data[i] = (local_data[i] + local_data[i+1] + local_data[i+2])/3;
		}

		location -= 2;
	}

	for(i=0;i<array_size;i++)
	{
		total += local_data[i];
	}

	MPI_Reduce(&total,&checksum,1,MPI_DOUBLE,MPI_SUM,0,comm);

	if (rank == 0)
	{
		time = MPI_Wtime() - start_time;
		checksum = checksum/(size*array_size);
		printf("Checksum = %.10lf\n",checksum);
		printf("Time taken = %lf seconds\n",time);
	}

	MPI_Finalize();

	return;
}
