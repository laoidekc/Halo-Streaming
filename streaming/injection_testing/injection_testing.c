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
	MPI_Status status;
	MPI_Request request;
	FILE *f;

	if(size != 2 && rank == 0)
	{
		printf("Must be run on exactly two processes\n");
	}

	if(argc < 4 && rank == 0)
	{
		printf("Need to supply number of messages, message size, and amount of calculation as arguments\n");
		return 0;
	}

	int num_messages = atoi(argv[1]);
	int message_size = atoi(argv[2]);
	int calculation_size = atoi(argv[3]);

	int i, flag, count;

	double send_times[num_messages];
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

	MPI_Request send_test_requests[num_messages];
	MPI_Status send_test_statuses[num_messages];

	// Send tests
	if(rank==1)
	{
		printf("Proc 1 cleaning up send tests\n");
		for(i=0;i<num_messages;i++)
		{
			MPI_Irecv(message_data,message_size,MPI_DOUBLE,0,0,comm,&send_test_requests[i]);
		}

		printf("Send test clean up complete\n");
	}

	
	MPI_Barrier(comm);

	if(rank==1)
	{
		MPI_Waitall(num_messages,send_test_requests,send_test_statuses);
	}

	if(rank==0)
	{
		char filename[20] = "sends_output.dat";
		f = fopen(filename, "w");
		send_times[0]=MPI_Wtime();
		for(i=0;i<num_messages;i++)
		{
			MPI_Send(message_data,message_size,MPI_DOUBLE,1,0,comm);
			calculation_delay(calculation_data,calculation_size);
			send_times[i+1]=MPI_Wtime();
		}

		fprintf(f,"Messages:\t%i\nMessage Size\t%i\nCalculation\t%i\nWalltime\t\tGap\n%lf\n",num_messages,message_size,calculation_size,send_times[0]);

		for(i=1;i<num_messages+1;i++)
		{
			fprintf(f,"%lf\t\t%lf\n",send_times[i],send_times[i]-send_times[i-1]);
		}

		fclose(f);
		printf("Send tests finished\n");
	}

	
	// Receive tests

	MPI_Barrier(comm);

	if(rank==0)
	{
		printf("Proc 0 setting up receive tests\n");
		count = 0;
		send_times[1] = send_times[0];
		while(send_times[1]-send_times[0]<10 && count<num_messages)
		{
			MPI_Isend(message_data,message_size,MPI_DOUBLE,1,0,comm,&request);
			send_times[0] = MPI_Wtime();
			send_times[1] = send_times[0];
			flag = 0;
			while(flag==0 && send_times[1]-send_times[0]<10)
			{
				MPI_Test(&request,&flag,&status);
				send_times[1] = MPI_Wtime();
			}
			count++;
			//printf("Message %i send after %lf seconds\n",count,send_times[1]-send_times[0]);
		}

		printf("Setup completed with %i messages\n", count);
	}

	MPI_Barrier(comm);

	if(rank==1)
	{	
		printf("Beginning receive tests\n");
		double receive_times[num_messages+1];
		char receive_filename[20] = "receives_output.dat";
		f = fopen(receive_filename, "w");
		receive_times[0] = MPI_Wtime();
		for(i=0;i<num_messages;i++)
		{
			MPI_Recv(message_data,message_size,MPI_DOUBLE,0,0,comm,&status);
			calculation_delay(calculation_data,calculation_size);
			receive_times[i+1] = MPI_Wtime();
		}
		fprintf(f,"Messages:\t%i\nMessage Size\t%i\nCalculation\t%i\nWalltime\t\tGap\n%lf\n",num_messages,message_size,calculation_size,receive_times[0]);

		for(i=1;i<num_messages+1;i++)
		{
			fprintf(f,"%lf\t\t%lf\n",receive_times[i],receive_times[i]-receive_times[i-1]);
		}

		fclose(f);
		printf("Receive tests finished\n");
	}

	MPI_Finalize();
}

void calculation_delay(double *calculation_data, int calculation_size)
{
	for(int i=0;i<calculation_size;i++)
	{
		// One calculation on each element of data
		calculation_data[i] = calculation_data[i] + 1;

		// Multiple calculations on one data element
		//calculation_data[0]++;
	}

}
