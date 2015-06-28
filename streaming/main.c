#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>

void initialise_data(double*, int);
void IO(int, int, int, double*, char*);

int main(int argc, char *argv[])
{
	MPI_Init(NULL, NULL);

	// Parameter declarations
	int array_size = 100;
	int iterations = 1000000;
	int send_buffer_size = 10;
	int receive_buffer_size = 10;
	int message_size = 50;
	int halo_size = 2;
	int num_runs = 10;

	// Output files
	char stream_filename[20] = "out.bin";
	char halo_filename[20] = "halo_out.bin";
	FILE *f;

	int rank, size, i, j, k, l, location_send, location_receive, send_tag, receive_tag, leftover, flag, send_iterations, receive_iterations, index;
	double *local_data, *new, *temporary, *send_buffer, *receive_buffer, start_time, time;
	size_t halo_bytes = halo_size*sizeof(double);

	local_data = malloc((array_size+halo_size)*sizeof(double));
	new = malloc((array_size+halo_size)*sizeof(double));
	send_buffer = malloc((halo_size*message_size*send_buffer_size)*sizeof(double));
	receive_buffer = malloc((halo_size*message_size*receive_buffer_size)*sizeof(double));

	MPI_Comm stream_comm = MPI_COMM_WORLD;
	MPI_Comm exchange_comm;
	MPI_Comm_dup(stream_comm,&exchange_comm);
	MPI_Comm_rank(stream_comm, &rank);
	MPI_Comm_size(stream_comm, &size);
	MPI_Request request_send[send_buffer_size], request_receive[receive_buffer_size], send_request_left, send_request_right, receive_request_left, receive_request_right;
	MPI_Status status_send[send_buffer_size], status_receive[receive_buffer_size], send_status_left, send_status_right, receive_status_left, receive_status_right;

	if(rank == 0)
	{
		f = fopen("data.txt", "w");
		fprintf(f,"Array size: %i\t\tIterations: %i\t\tMessage size: %i\t\tProcessors used: %i\n",array_size,iterations,message_size,size);
		fprintf(f,"Stream time\tExchange time\n");
	}

	for(l=0;l<num_runs;l++)
	{
		location_send = array_size - halo_size;
		location_receive = array_size - halo_size;
		send_iterations = 0;
		receive_iterations = 0;

		// Number of iterations that can be completed without communication
		int triangle_iterations = (int)fmin((double)array_size/halo_size,(double)iterations);

		// Find neighbours
		int neighbour_left = (rank + size - 1) % size;
		int neighbour_right = (rank + 1) % size;

		// Position in global array
		int start_point = rank*array_size;

		// Initialise local data
		initialise_data(local_data, array_size);

		// Initialise persistent sends
		for(send_tag=0;send_tag<send_buffer_size;send_tag++)
		{
			MPI_Send_init(&send_buffer[halo_size*message_size*send_tag],halo_size*message_size,MPI_DOUBLE,neighbour_left,send_tag,stream_comm,&request_send[send_tag]);	
		}

		// Initialise non-blocking receives
		for(receive_tag=0;receive_tag<receive_buffer_size;receive_tag++)
		{
			MPI_Recv_init(&receive_buffer[halo_size*message_size*receive_tag],halo_size*message_size,MPI_DOUBLE,neighbour_right,receive_tag,stream_comm,&request_receive[receive_tag]);
			MPI_Start(&request_receive[receive_tag]);
		}

		// Reset tags
		send_tag = 0;
		receive_tag = 0;

		MPI_Barrier(stream_comm);

		// Begin halo-streaming timing
		if(rank == 0)
		{
			start_time = MPI_Wtime();
		}

		// Compute initial triangle
		while(send_iterations<triangle_iterations)
		{
			// Wait for request to become available
			MPI_Wait(&request_send[send_tag],&status_send[send_tag]);

			// Loop over all data to be sent	
			for(j=0;j<message_size;j++)
			{
				// Copy data to send buffer 
				memcpy(&send_buffer[halo_size*message_size*send_tag+j*halo_size],local_data,halo_bytes);

				// Update all possible values (in place)
				for(i=0;i<location_send;i++)
				{
					local_data[i] = (local_data[i] + local_data[i+1] + local_data[i+2])/3;
				}

				// Reduce the number of values that can be updated during the next iteration
				location_send -= halo_size;
			}

			// Start new send
			MPI_Start(&request_send[send_tag]);
			
			// Change tag
			send_tag = (send_tag + 1) % send_buffer_size;
			send_iterations += message_size;

			// Check if receive has completed
			MPI_Test(&request_receive[receive_tag],&flag,&status_receive[receive_tag]);
		
			if(flag)
			{
				// Loop over all received data
				for(j=0;j<message_size;j++)
				{
					// Copy data from receive buffer
					memcpy(&local_data[array_size],&receive_buffer[halo_size*message_size*receive_tag+j*halo_size],halo_bytes);	

					// Update values on the edge of the triangle
					for(i=0;i+receive_iterations<send_iterations;i++)
					{
						for(k=0;k<halo_size;k++)
						{
							index = array_size + k - halo_size*(i+1);
							local_data[index] = (local_data[index] + local_data[index+1] + local_data[index+2])/3;
						}
					}
					receive_iterations++;
				}

				// Increase size of triangle that can be computed
				location_send += halo_size*message_size;
				triangle_iterations += message_size;
				if(triangle_iterations>iterations)
				{
					triangle_iterations = iterations;
				}

				// Begin new receive
				MPI_Start(&request_receive[receive_tag]);

				// Change tag
				receive_tag = (receive_tag + 1) % receive_buffer_size;
			}
		}

		// Extend triangle until its peak reaches the maximum number of iterations
		while(send_iterations<iterations)
		{
			// Wait for requests
			MPI_Wait(&request_receive[receive_tag],&status_receive[receive_tag]);
			MPI_Wait(&request_send[send_tag],&status_send[send_tag]);

			// Loop over all received data
			for(j=0;j<message_size;j++)
			{
				// Copy from receive buffer
				memcpy(&local_data[array_size],&receive_buffer[halo_size*message_size*receive_tag+j*halo_size],halo_bytes);

				// Compute diagonal
				for(i=array_size-halo_size;i>=0;i-=halo_size)
				{
					for(k=0;k<halo_size;k++)
					{
						index = i+k;
						local_data[index] = (local_data[index] + local_data[index+1] + local_data[index+2])/3;
					}
				}

				// Copy to send buffer
				memcpy(&send_buffer[halo_size*message_size*send_tag+j*halo_size],local_data,halo_bytes);
			}

			// Begin new send and receive
			MPI_Start(&request_send[send_tag]);
			MPI_Start(&request_receive[receive_tag]);

			// Increase tags
			receive_tag = (receive_tag + 1) % receive_buffer_size;
			send_tag = (send_tag + 1) % send_buffer_size;

			receive_iterations += message_size;
			send_iterations += message_size;
		}

		// Fill in remaining inverted triangle
		while(receive_iterations<iterations)
		{
			// Wait for receive to complete
			MPI_Wait(&request_receive[receive_tag], &status_receive[receive_tag]);

			// Loop over all received data
			for(j=0;j<message_size;j++)
			{
				// Copy data from receive buffer
				memcpy(&local_data[array_size],&receive_buffer[halo_size*message_size*receive_tag+j*halo_size],halo_bytes);

				// Update all possible values (in place)
				for(i=location_receive;i<array_size;i++)
				{
					local_data[i] = (local_data[i] + local_data[i+1] + local_data[i+2])/3;
				}

				// Increase the number of values to be updated during the next iteration
				location_receive -= halo_size;
			}

			// Begin new receive
			MPI_Start(&request_receive[receive_tag]);

			// Change tag
			receive_tag = (receive_tag + 1) % receive_buffer_size;

			receive_iterations += message_size;
		}

		MPI_Barrier(stream_comm);

		// End halo-streaming timing
		if (rank == 0)
		{
			time = MPI_Wtime() - start_time;
			fprintf(f,"%f\t",time);
		}

		// Find new position in global array
		start_point = (start_point + iterations) % (size*array_size);

		// I/O operations
		IO(start_point, array_size, size, local_data, stream_filename);

		// Cancel all outstanding communications
		for(receive_tag=0;receive_tag<receive_buffer_size;receive_tag++)
		{
			MPI_Cancel(&request_receive[receive_tag]);
		}

		// Re-initialise data for halo-exchange
		initialise_data(&local_data[halo_size/2],array_size);

		MPI_Barrier(exchange_comm);

		// Begin timing for halo-exchange
		if(rank == 0)
		{
			start_time = MPI_Wtime();
		}

		for(j=0;j<iterations;j++)
		{
			// Send and receive halos
			MPI_Isend(&local_data[1],1,MPI_DOUBLE,neighbour_left,j,exchange_comm,&send_request_left);
			MPI_Isend(&local_data[array_size],1,MPI_DOUBLE,neighbour_right,j,exchange_comm,&send_request_right);

			MPI_Irecv(&local_data[0],1,MPI_DOUBLE,neighbour_left,j,exchange_comm,&receive_request_left);
			MPI_Irecv(&local_data[array_size+1],1,MPI_DOUBLE,neighbour_right,j,exchange_comm,&receive_request_right);

			// Update data
			for(i=2;i<=array_size-1;i++)
			{
				new[i] = (local_data[i-1] + local_data[i] + local_data[i+1])/3;
			}

			MPI_Wait(&receive_request_left,&receive_status_left);
			MPI_Wait(&receive_request_right,&receive_status_right);

			new[1] = (local_data[0] + local_data[1] + local_data[2])/3;
			new[array_size] = (local_data[array_size-1] + local_data[array_size] + local_data[array_size+1])/3;

			MPI_Wait(&send_request_left,&send_status_left);
			MPI_Wait(&send_request_right,&send_status_right);

			// Swapping pointers to arrays
			temporary = new;
			new = local_data;
			local_data = temporary;
		}

		MPI_Barrier(exchange_comm);

		// End halo-exchange timing
		if (rank == 0)
		{
			time = MPI_Wtime() - start_time;
			fprintf(f,"%f\n",time);
		}

		// Find position in global array
		start_point = array_size*rank;

		// I/O operations
		IO(start_point,array_size,size,&local_data[halo_size/2],halo_filename);
	}

	// Free arrays
	free(local_data);
	free(new);
	free(send_buffer);
	free(receive_buffer);

	if(rank == 0)
	{
		fclose(f);
	}

	MPI_Finalize();

	return 0;
}

void initialise_data(double *local_data, int array_size)
{
	srand(227);
	int i;

	for(i=0;i<array_size;i++)
	{
		local_data[i] = rand()%10000;
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
		// Rearrange data before writing to file
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

		// Createsplit datatype
		int blocklengths[2] = {leftover,array_size-leftover};
		int displacements[2] = {0,(size - 1)*array_size + leftover};
		MPI_Type_indexed(2,blocklengths,displacements,MPI_DOUBLE,&filetype);
	}

	else
	{
		disp = start_point*sizeof(double);

		// Create contiguous datatype
		MPI_Type_contiguous(array_size,MPI_DOUBLE,&filetype);
	}

	MPI_Type_commit(&filetype);

	MPI_File_open(comm,filename,amode,MPI_INFO_NULL,&fh);

	MPI_File_set_view(fh,disp,MPI_DOUBLE,filetype,"native",MPI_INFO_NULL);

	MPI_File_write_all(fh,local_data,array_size,MPI_DOUBLE,&file_status);

	MPI_File_close(&fh);

	MPI_Type_free(&filetype);
}
