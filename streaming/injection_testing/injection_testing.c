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
	char injection_filename[22] = "injection_output.dat";
	char extraction_filename[22] = "extraction_output.dat";
	FILE *f;
	int i;

	if(size != 2 && rank == 0)
	{
		printf("Must be run on exactly two processes\n");
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

	MPI_Request send_test_requests[num_messages];
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

	// Injection test
	if(rank==1)
	{
		printf("Proc 1 posting receives to prepare for injection test\n");
		for(i=0;i<num_messages;i++)
		{
			MPI_Irecv(message_data,message_size,MPI_DOUBLE,0,0,comm,&send_test_requests[i]);
		}
		printf("Proc 1 preparation for injection test complete\n");
	}

	MPI_Barrier(comm);

	if(rank==1)
	{
		MPI_Waitall(num_messages,send_test_requests,MPI_STATUSES_IGNORE);
		printf("Proc 1 injection test complete\n");
	}

	if(rank==0)
	{
		printf("Proc 0 beginning injection test\n");
		send_times[0]=MPI_Wtime();
		for(i=0;i<num_messages;i++)
		{
			MPI_Send(message_data,message_size,MPI_DOUBLE,1,0,comm);
			calculation_delay(calculation_data,calculation_size);
			send_times[i+1]=MPI_Wtime();
		}

		f = fopen(injection_filename, "w");
		fprintf(f,"#Messages:\t%i\n#Message Size\t%i\n#Calculation\t%i\n#Index\t\tWalltime\t\tGap\n#0\t\t%lf\n",num_messages,message_size,calculation_size,send_times[0]);

		for(i=1;i<num_messages+1;i++)
		{
			fprintf(f,"%i\t\t%lf\t\t%lf\n",i,send_times[i],send_times[i]-send_times[i-1]);
		}
		fclose(f);

		printf("Proc 0 injection test complete\n");
	}

	
	// Receive tests

	MPI_Barrier(comm);

	if(rank==0)
	{
		printf("Proc 0 posting sends to prepare for extraction test\n");
		for(i=0;i<num_messages;i++)
		{
			MPI_Isend(message_data,message_size,MPI_DOUBLE,1,0,comm,&send_test_requests[i]);
		}
		printf("Proc 0 preparation for extraction test complete\n");
	}

	MPI_Barrier(comm);

	if(rank==0)
	{
		MPI_Waitall(num_messages,send_test_requests,MPI_STATUSES_IGNORE);
		printf("Proc 0 extraction test complete\n");
	}

	if(rank==1)
	{	
		printf("Proc 1 beginning extraction test\n");
		double receive_times[num_messages+1];
		receive_times[0] = MPI_Wtime();
		for(i=0;i<num_messages;i++)
		{
			MPI_Recv(message_data,message_size,MPI_DOUBLE,0,0,comm,MPI_STATUS_IGNORE);
			calculation_delay(calculation_data,calculation_size);
			receive_times[i+1] = MPI_Wtime();
		}

		f = fopen(extraction_filename, "w");
		fprintf(f,"#Messages:\t%i\n#Message Size\t%i\n#Calculation\t%i\n#Index\t\tWalltime\t\tGap\n#0\t\t%lf\n",num_messages,message_size,calculation_size,receive_times[0]);

		for(i=1;i<num_messages+1;i++)
		{
			fprintf(f,"%i\t\t%lf\t\t%lf\n",i,receive_times[i],receive_times[i]-receive_times[i-1]);
		}
		fclose(f);

		printf("Proc 1 extraction test complete\n");
	}

	MPI_Finalize();
}

void calculation_delay(double *calculation_data, int calculation_size)
{
	for(int i=0;i<calculation_size;i++)
	{
		// One calculation on each element of data
		//calculation_data[i] = calculation_data[i] + 1;

		// Multiple calculations on one data element
		calculation_data[0]++;
	}

}
