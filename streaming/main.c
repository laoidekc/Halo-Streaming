#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <stdlib.h>

void initialise_data(double*, int);
void IO(int, int, int, double*, char*);

int main(int argc, char *argv[])
{
	MPI_Init(NULL, NULL);


	// Parameter declarations
	int array_size = 100000;
	int iterations = 10000;
	int send_buffer_size = 10;
	int message_size = 10;
	int halo_size = 2;

	// Output files
	char stream_filename[20] = "out.bin";
	char halo_filename[20] = "halo_out.bin";

	int rank, size, i, j, location_send = array_size - halo_size, location_receive = array_size - halo_size, send_tag, receive_tag, leftover, flag;
	double *local_data, *new, *temporary, *send_buffer, *receive_buffer, start_time, time;

	int send_buffer_space = send_buffer_size;
	int receive_buffer_space = send_buffer_size;
	int send_iterations = 0;
	int receive_iterations = 0;

	local_data = malloc((array_size+halo_size)*sizeof(double));
	new = malloc((array_size+halo_size)*sizeof(double));
	send_buffer = malloc((halo_size*message_size*send_buffer_size)*sizeof(double));
	receive_buffer = malloc((halo_size*message_size*send_buffer_size)*sizeof(double));

	MPI_Comm stream_comm = MPI_COMM_WORLD;
	MPI_Comm exchange_comm;
	MPI_Comm_dup(stream_comm,&exchange_comm);
	MPI_Comm_rank(stream_comm, &rank);
	MPI_Comm_size(stream_comm, &size);
	MPI_Request request_send[send_buffer_size], request_receive[send_buffer_size], request_left, request_right;
	MPI_Status status_send[send_buffer_size], status_receive[send_buffer_size], status_left, status_right;

	// Finding neighbours
	int neighbour_left = (rank + size - 1) % size;
	int neighbour_right = (rank + 1) % size;

	// Position in global array
	int start_point = rank*array_size;

	// Initialising local data
	initialise_data(local_data, array_size);

	// Initialising persistent sends
	for(send_tag=0;send_tag<send_buffer_size;send_tag++)
	{
		MPI_Send_init(&send_buffer[halo_size*message_size*send_tag],halo_size*message_size,MPI_DOUBLE,neighbour_left,send_tag,stream_comm,&request_send[send_tag]);	
	}

	// Initialising non-blocking receives
	for(receive_tag=0;receive_tag<send_buffer_size;receive_tag++)
	{
		MPI_Irecv(&receive_buffer[halo_size*message_size*receive_tag],halo_size*message_size,MPI_DOUBLE,neighbour_right,receive_tag,stream_comm,&request_receive[receive_tag]);
	}

	send_tag = 0;
	receive_tag = 0;

	MPI_Barrier(stream_comm);

	// Beginning halo-streaming timing
	if(rank == 0)
	{
		start_time = MPI_Wtime();
	}

	for(send_iterations=0;send_iterations<iterations;send_iterations+=message_size)
	{
		// Waiting for request to become available
		MPI_Wait(&request_send[send_tag],&status_send[send_tag]);
		
		for(j=0;j<message_size;j++)
		{
			// Copying data to send buffer 
			for(i=0;i<halo_size;i++)
			{
				send_buffer[halo_size*message_size*send_tag+j*halo_size+i] = local_data[i];
			}

			// Updating all possible values (in place)
			for(i=0;i<location_send;i++)
			{
				local_data[i] = (local_data[i] + local_data[i+1] + local_data[i+2])/3;
			}

			// Reducing the number of values that can be updated during the next iteration
			location_send -= halo_size;
		}

		MPI_Start(&request_send[send_tag]);
			
		// Changing tag
		send_tag = (send_tag + 1) % send_buffer_size;

		MPI_Test(&request_receive[receive_tag],&flag,&status_receive[receive_tag]);
		
		if(flag)
		{
			for(j=0;j<message_size;j++)
			{
				for(i=0;i<halo_size;i++)
				{
					local_data[array_size+i] = receive_buffer[halo_size*message_size*receive_tag+j*halo_size+i];
				}

				for(i=location_receive;i<array_size;i++)
				{
					local_data[i] = (local_data[i] + local_data[i+1] + local_data[i+2])/3;
				}

				location_receive -= halo_size;
			}

			MPI_Irecv(&receive_buffer[halo_size*message_size*receive_tag],halo_size*message_size,MPI_DOUBLE,neighbour_right,receive_tag,stream_comm,&request_receive[receive_tag]);

			receive_tag = (receive_tag + 1) % send_buffer_size;

			receive_iterations += message_size;
		}
	}
	
	// Waiting for all sends to complete (this may be unnecessary)
	MPI_Waitall(send_buffer_size,request_send,status_send);

	// Finding new position in global array
	start_point = (start_point + iterations) % (size*array_size);

	for(;receive_iterations<iterations;receive_iterations+=message_size)
	{

		// Waiting for receive to complete
		MPI_Wait(&request_receive[receive_tag], &status_receive[receive_tag]);

		for(j=0;j<message_size;j++)
		{
			for(i=0;i<halo_size;i++)
			{
				local_data[array_size+i] = receive_buffer[halo_size*message_size*receive_tag+j*halo_size+i];
			}

			for(i=location_receive;i<array_size;i++)
			{
				local_data[i] = (local_data[i] + local_data[i+1] + local_data[i+2])/3;
			}

			location_receive -= halo_size;
		}

		MPI_Irecv(&receive_buffer[halo_size*message_size*receive_tag],halo_size*message_size,MPI_DOUBLE,neighbour_right,receive_tag,stream_comm,&request_receive[receive_tag]);

		receive_tag = (receive_tag + 1) % send_buffer_size;
	}

	MPI_Barrier(stream_comm);

	// Ending halo-streaming timing
	if (rank == 0)
	{
		time = MPI_Wtime() - start_time;
		printf("Halo streaming time taken = %lf seconds\n",time);
	}

	IO(start_point, array_size, size, local_data, stream_filename);

	// Re-initialising data for halo-exchange
	initialise_data(&local_data[halo_size/2],array_size);

	MPI_Barrier(exchange_comm);

	// Beginning timing for halo-exchange
	if(rank == 0)
	{
		start_time = MPI_Wtime();
	}

	for(j=0;j<iterations;j++)
	{
		// Sending and receiving halos
		MPI_Issend(&local_data[1],1,MPI_DOUBLE,neighbour_left,j,exchange_comm,&request_left);
		MPI_Issend(&local_data[array_size],1,MPI_DOUBLE,neighbour_right,j,exchange_comm,&request_right);

		MPI_Recv(&local_data[array_size+1],1,MPI_DOUBLE,neighbour_right,j,exchange_comm,&status_right);
		MPI_Recv(&local_data[0],1,MPI_DOUBLE,neighbour_left,j,exchange_comm,&status_left);

		MPI_Wait(&request_left,&status_left);
		MPI_Wait(&request_right,&status_right);

		// Updating data
		for(i=1;i<=array_size;i++)
		{
			new[i] = (local_data[i-1] + local_data[i] + local_data[i+1])/3;
		}

		// Swapping pointers to arrays
		temporary = new;
		new = local_data;
		local_data = temporary;
	}

	MPI_Barrier(exchange_comm);

	// Ending halo-exchange timing
	if (rank == 0)
	{
		time = MPI_Wtime() - start_time;
		printf("Halo exchange time taken = %lf seconds\n",time);
	}

	start_point = array_size*rank;

	IO(start_point,array_size,size,&local_data[halo_size/2],halo_filename);

	// Freeing arrays
	free(local_data);
	free(new);
	free(send_buffer);
	free(receive_buffer);

	MPI_Finalize();

	return;
}

void initialise_data(double *local_data, int array_size)
{
	int i;

	for(i=0;i<array_size;i++)
	{
		local_data[i] = (double)i;
	}
}

void IO(int start_point, int array_size, int size, double *local_data, char *filename)
{
	int amode = MPI_MODE_CREATE | MPI_MODE_WRONLY;
	MPI_File fh;
	MPI_Datatype filetype;
	MPI_Offset disp;
	MPI_Comm comm = MPI_COMM_WORLD;
	int leftover = start_point + array_size - array_size*size;
	int i;
	MPI_Status file_status;

	
	if(leftover > 0)
	{
		// Rearranging data before writing to file
		disp = 0;
		double temp[leftover];
		for(i=0;i<leftover;i++)
		{
			temp[i] = local_data[array_size-leftover+i];
		}

		for(i=array_size-1;i>=leftover;i--)
		{
			local_data[i] = local_data[i - leftover];
		}

		for(i=0;i<leftover;i++)
		{
			local_data[i] = temp[i];
		}

		// Creating split datatype
		int blocklengths[2] = {leftover,array_size-leftover};
		int displacements[2] = {0,(size - 1)*array_size + leftover};
		MPI_Type_indexed(2,blocklengths,displacements,MPI_DOUBLE,&filetype);
	}

	else
	{
		disp = start_point*sizeof(double);

		// Creating contiguous datatype
		MPI_Type_contiguous(array_size,MPI_DOUBLE,&filetype);
	}

	MPI_Type_commit(&filetype);

	MPI_File_open(comm,filename,amode,MPI_INFO_NULL,&fh);

	MPI_File_set_view(fh,disp,MPI_DOUBLE,filetype,"native",MPI_INFO_NULL);

	MPI_File_write_all(fh,local_data,array_size,MPI_DOUBLE,&file_status);

	MPI_File_close(&fh);

	MPI_Type_free(&filetype);
}
