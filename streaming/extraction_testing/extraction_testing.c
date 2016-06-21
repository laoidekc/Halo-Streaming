#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

void calculation_delay(double*, int);

//This code is for testing the extraction and data retrieval rates of MPI
int main(int argc, char *argv[])
{
	MPI_Init(NULL,NULL);
	
	int rank, size;
	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &size);
	char extraction_filename[128] = "extraction_output";
	time_t current_time;
	struct tm * time_info;
	char time_string[21];
	time(&current_time);
	time_info = localtime(&current_time);
	strftime(time_string,sizeof(time_string),"_%F_%H-%M-%S",time_info);

	strcat(extraction_filename,"_NM");
	strcat(extraction_filename,argv[1]);
	strcat(extraction_filename,"_MS");
	strcat(extraction_filename,argv[2]);
	strcat(extraction_filename,"_CS");
	strcat(extraction_filename,argv[3]);
	strcat(extraction_filename,time_string);
	strcat(extraction_filename,".dat");

	FILE *f;
	int i;

	if(size < 2 && rank == 0)
	{
		printf("Must be run on two or more processes\n");
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

	// Receive tests

	MPI_Barrier(comm);

	if(rank==0)
	{
		printf("Proc 0 posting sends to prepare for extraction test\n");
		for(i=0;i<num_messages;i++)
		{
			MPI_Isend(message_data,message_size,MPI_DOUBLE,1,0,comm,&requests[i]);
		}
		printf("Proc 0 preparation for extraction test complete\n");
	}

	MPI_Barrier(comm);

	if(rank==0)
	{
		MPI_Waitall(num_messages,requests,MPI_STATUSES_IGNORE);
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
