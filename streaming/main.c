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
	int array_size = 983040000;		// Number of doubles per processor
	int iterations = 10000;			// Length of simulation
	int send_buffer_size = 1000;	// Number of requests that can be used to send data. Should be greater than the maximum number of expected outstanding messages
	int receive_buffer_size = 1000;	// Same thing for the receive buffer
	int message_size = 50;			// Number of iterations that are completed before performing communications. This number should divide evenly into both array_size and iterations. It also must be less than array_size/2.
	int halo_size = 2;				// Number of data points contained in a processor's halo region. This is the sum of the halo regions in both directions.
	int num_runs = 10;				// Number of trials the program will perform for both halo streaming and halo exchange.
	int halo_depth = atoi(argv[1]);                            // Number of data points the halo_exchange code will transfer. This should be 1 for normal halo exchange.

	// Output files
	char stream_filename[20] = "out.bin";		// Binary output of the final halo streaming data
	char halo_filename[20] = "halo_out.bin";	// Binary output of the final halo exchange data
	char output_filename[20] = "data";			// Will contain the timing results
       	strcat(output_filename, argv[1]);			
	strcat(output_filename, ".txt");			

	FILE *f;

	int rank, size, i, j, k, l, location_send, location_receive, send_tag, receive_tag, leftover, flag, send_iterations, receive_iterations, index;
	double *local_data, *exchange_array, *new, *temporary, *send_buffer, *receive_buffer, start_time, stream_time[num_runs], exchange_time[num_runs], average_time = 0, deviation_time = 0;

	size_t halo_bytes = halo_size*sizeof(double);	// Size of halo region. Will be used later for various memcpys.

	MPI_Comm stream_comm = MPI_COMM_WORLD;		// Communicator for all streaming communications
	MPI_Comm exchange_comm;						// Communicator for all exchange communications
	MPI_Comm_dup(stream_comm,&exchange_comm);

	MPI_Comm_rank(stream_comm, &rank);		// Assigns rank to each processor
	MPI_Comm_size(stream_comm, &size);		// Number of processors involved
	MPI_Request request_send[send_buffer_size], request_receive[receive_buffer_size], send_request_left, send_request_right, receive_request_left, receive_request_right;	// Request arrays are used for streaming, the rest for halo exchange.
	MPI_Status status_send[send_buffer_size], status_receive[receive_buffer_size], send_status_left, send_status_right, receive_status_left, receive_status_right;		// Status arrays are used for streaming, the rest for halo exchange.

       	array_size = array_size/size;

	local_data = malloc((array_size+halo_size)*sizeof(double));		// Primary work buffer for halo streaming. Has enough space for the starting data plus one halo region.
	exchange_array = malloc((array_size+halo_size*halo_depth)*sizeof(double)); // Primary work buffer for halo exchange.
	new = malloc((array_size+halo_size*halo_depth)*sizeof(double));	// Secondary work buffer to be used in halo exchange.
	send_buffer = malloc((halo_size*message_size*send_buffer_size)*sizeof(double));		// Send buffer for halo streaming
	receive_buffer = malloc((halo_size*message_size*receive_buffer_size)*sizeof(double));	// Receive buffer for halo streaming


	// Master processor writes parameters to timing file
	if(rank == 0)
	{
		f = fopen(output_filename, "w");
		fprintf(f,"Array size: %i\t\tIterations: %i\t\tMessage size: %i\t\tHalo depth: %i\tProcessors used: %i\n",array_size,iterations,message_size,halo_depth,size);
		fprintf(f,"Run\t\tStream time\tExchange time\n");
	}

	// Loop over number of trials
	for(l=0;l<num_runs;l++)
	{
		location_send = array_size - halo_size;		// Denotes the point in the array past which the processor no longer has enough information to perform updates. Initialised to array_size - halo_size because there will be halo_size number of points that depend on data from the neighbouring processor.
		location_receive = array_size - halo_size; // Denotes the start of the area that can be calculated when new data arrives.

		send_iterations = 0;	// Number of iterations worth of data that have been sent
		receive_iterations = 0;	// Number of iterations worth of data that have been received

		// Number of iterations that can be completed without communication
		int triangle_iterations = (int)fmin((double)array_size/halo_size,(double)iterations);

		// Find neighbours
		int neighbour_left = (rank + size - 1) % size;
		int neighbour_right = (rank + 1) % size;

		// Position in global array. This tracks how far the domain of a given processor shifts due to halo streaming.
		int start_point = rank*array_size;

		// Initialise local data
		initialise_data(local_data, array_size);

		// Initialise persistent sends
		for(send_tag=0;send_tag<send_buffer_size;send_tag++)
		{
			MPI_Send_init(&send_buffer[halo_size*message_size*send_tag],halo_size*message_size,MPI_DOUBLE,neighbour_left,send_tag,stream_comm,&request_send[send_tag]);	
		}

		// Initialise persistent receives
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
		/*	while(send_iterations<triangle_iterations)
		{
			// Wait for request to become available
			MPI_Wait(&request_send[send_tag],&status_send[send_tag]);

			// Completes iterations until a block of data to be sent is complete	
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

			// Updates the number of iterations that have been sent
			send_iterations += message_size;

			// Check if receive has completed
			MPI_Test(&request_receive[receive_tag],&flag,&status_receive[receive_tag]);
		
			if(flag)
			{
				// Contines until all data from the message has been processed
				for(j=0;j<message_size;j++)
				{
					// Copy data from receive buffer
					memcpy(&local_data[array_size],&receive_buffer[halo_size*message_size*receive_tag+j*halo_size],halo_bytes);	

					// Update all values that depend on the new data. Stops when it catches up with the data that has been processed for sending
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

				// Increase size of triangle that can be computed, due to new information having arrived
				location_send += halo_size*message_size;
				triangle_iterations += message_size;

				// Ensures that processor does not calculate past the final iteration
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

				// Calculate all values that depend on the new data
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

			// Update iterations
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

			// Update iterations
			receive_iterations += message_size;
			}*/

		MPI_Barrier(stream_comm);

		// End halo-streaming timing and write time to file
		if (rank == 0)
		{
			stream_time[l] = MPI_Wtime() - start_time;
			fprintf(f,"%i\t\t%f\t",l,stream_time[l]);
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
		initialise_data(&exchange_array[halo_size*halo_depth/2],array_size);

		MPI_Barrier(exchange_comm);

		// Begin timing for halo-exchange
		if(rank == 0)
		{
			start_time = MPI_Wtime();
		}

	    for(j=0;j<iterations;j+=halo_depth)
		{
			// Send and receive halos
			MPI_Isend(&exchange_array[halo_depth],halo_depth,MPI_DOUBLE,neighbour_left,j,exchange_comm,&send_request_left);
			MPI_Isend(&exchange_array[array_size],halo_depth,MPI_DOUBLE,neighbour_right,j,exchange_comm,&send_request_right);

			MPI_Irecv(&exchange_array[0],halo_depth,MPI_DOUBLE,neighbour_left,j,exchange_comm,&receive_request_left);
			MPI_Irecv(&exchange_array[array_size+halo_depth],halo_depth,MPI_DOUBLE,neighbour_right,j,exchange_comm,&receive_request_right);

			// Wait until data has been received
			MPI_Wait(&receive_request_left,&receive_status_left);
			MPI_Wait(&receive_request_right,&receive_status_right);
			// Wait until data has been sent
			MPI_Wait(&send_request_left,&send_status_left);
			MPI_Wait(&send_request_right,&send_status_right);

			for(k=1;k<=halo_depth;k++)
			{  
			  // Update central data
			  for(i=k;i<array_size+halo_size*halo_depth-k;i++)
			    {
			      new[i] = (exchange_array[i-1] + exchange_array[i] + exchange_array[i+1])/3;
			    }

			  // Update edge data
			  // new[1] = (local_data[0] + local_data[1] + local_data[2])/3;
			  // new[array_size] = (local_data[array_size-1] + local_data[array_size] + local_data[array_size+1])/3;

			  // Swapping pointers to arrays
			  temporary = new;
			  new = exchange_array;
			  exchange_array = temporary;
			}
		}

		MPI_Barrier(exchange_comm);

		// End halo-exchange timing and write time to file
		if(rank == 0)
		{
			exchange_time[l] = MPI_Wtime() - start_time;
			fprintf(f,"%f\n",exchange_time[l]);
		}

		// Find position in global array
		start_point = array_size*rank;

		// I/O operations
		IO(start_point,array_size,size,&exchange_array[halo_size*halo_depth/2],halo_filename);
	}

	if(rank == 0)
	{
		for(l=0;l<num_runs;l++)
		{
			average_time += stream_time[l];
		}
		average_time = average_time/num_runs;

		for(l=0;l<num_runs;l++)
		{
			deviation_time += (stream_time[l]-average_time)*(stream_time[l]-average_time);
		}

		deviation_time = sqrt(deviation_time/num_runs);

		fprintf(f,"Average:\t%f\t",average_time);

		average_time = 0;

		for(l=0;l<num_runs;l++)
		{
			average_time += exchange_time[l];
		}
		average_time = average_time/num_runs;

		fprintf(f,"%f\n",average_time);
		fprintf(f,"Deviation:\t%f\t",deviation_time);
		deviation_time = 0;

		for(l=0;l<num_runs;l++)
		{
			deviation_time += (exchange_time[l]-average_time)*(exchange_time[l]-average_time);
		}

		deviation_time = sqrt(deviation_time/num_runs);

		fprintf(f,"%f\n",deviation_time);
	}

	// Free arrays
	free(local_data);
	free(exchange_array);
	free(new);
	free(send_buffer);
	free(receive_buffer);

	// Close timing file
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
