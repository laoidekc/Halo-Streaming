#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int main(int argc, char *argv[])
{
	MPI_Init(NULL,NULL);
	
	int rank, size;
	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &size);
	char output_filename[22] = "basic_streams_output.dat";
	FILE *f;
	int i;

	if(size != 3 && rank == 0)
	{
		printf("Must be run on exactly three processes\n");
		return 0;
	}

	if(argc < 4 && rank == 0)
	{
		printf("Need to supply number of messages, message size, and amount of calculation as arguments\n");
		return 0;
	}

	int num_messages = atoi(argv[1]);
	int message_size = atoi(argv[2]);
	int calculation_size = atoi(argv[3]);

	MPI_Request requests[num_messages];
	double times[num_messages+1];

	double message_data[message_size];
	for(i=0;i<message_size;i++)
	{
		message_data[i] = 0;
	}

	double calculation_data[calculation_size];
	for(i=0;i<calculation_size;i++)
	{
		calculation_data[i] = rand()%10000;
	}

	if(rank==0)
	{
		printf("Proc 0 preparing receives\n");
		for(i=0;i<num_messages;i++)
		{
			MPI_Irecv(message_data,message_size,MPI_DOUBLE,1,0,comm,&requests[i]);
		}
		printf("Proc 0 preparation complete\n");
	}

	if(rank==2)
	{
		printf("Proc 2 preparing sends\n");
		for(i=0;i<num_messages;i++)
		{
			MPI_Isend(message_data,message_size,MPI_DOUBLE,1,0,comm,&requests[i]);

		}
		printf("Proc 2 preparation complete\n");
	}

	MPI_Barrier(comm);

	if(rank==1)
	{
		printf("Proc 1 beginning tests\n");
		times[0] = MPI_Wtime();
		for(i=0;i<num_messages;i++)
		{
			MPI_Recv(message_data,message_size,MPI_DOUBLE,2,0,comm,MPI_STATUS_IGNORE);
			// calculation here
			MPI_Send(message_data,message_size,MPI_DOUBLE,0,0,comm);
			times[i+1] = MPI_Wtime();
		}
		printf("Proc 1 finished tests and beginning I/O\n");

		f = fopen(output_filename, "w");
		fprintf(f,"#Messages:\t%i\n#Message Size\t%i\n#Calculation\t%i\n#Index\t\tWalltime\t\tGap\n#0\t\t%lf\n",num_messages,message_size,calculation_size,times[0]);
		for(i=1;i<num_messages+1;i++)
		{
			fprintf(f,"%i\t\t%lf\t\t%lf\n",i,receive_times[i],times[i]-times[i-1]);
		}
		fclose(f);
		printf("Proc 1 finished I/O"\n);
	}

	else
	{
		printf("Proc %i beginning cleanup\n",rank);
		MPI_Waitall(num_messages,requests,MPI_STATUSES_IGNORE);
		printf("Proc %i finished cleanup\n",rank);
	}



	MPI_Finalize();
}

