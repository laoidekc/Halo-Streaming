#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>

void initialise_data(double*, int);
void IO(int, int, int, double*, char*);
void outstanding_receives(int, MPI_Request*, int*, int*, int*, MPI_Status*, int*);

int main(int argc, char *argv[])
{
	MPI_Init(NULL, NULL);

	// Parameter declarations
	int array_size = 10000;			// Number of doubles per processor
	int per_process = 1;			// Set equal to 1 if array_size is the number of data elements per process. Set equal to 0 if array_size is equal to the global number of elements
	int iterations = 1000;			// Length of simulation
	int send_buffer_size = 100;	// Number of requests that can be used to send data. Should be greater than the maximum number of expected outstanding messages
	int receive_buffer_size = 100;	// Same thing for the receive buffer
	int message_size = 50;			// Number of iterations that are completed before performing communications. This number should divide evenly into both array_size and iterations. It also must be less than array_size/2.
	int num_runs = 0;				// Number of trials the program will perform for both halo streaming and halo exchange.

	int halo_size = 2;				// Number of data points contained in a processor's halo region. This is the sum of the halo regions in both directions.

	// Output files
	char stream_filename[20] = "streaming_out.bin";		// Binary output of the final halo streaming data
	char halo_filename[20] = "exchange_out.bin";	// Binary output of the final halo exchange data
	char output_filename[20] = "data";			// Will contain the timing results
	strcat(output_filename, argv[1]);			// File identifier
	strcat(output_filename, ".txt");			// File extension

	FILE *f;

	int rank, size, i, j, k, l, location_send, location_receive, send_tag, receive_tag, leftover, flag, send_iterations, receive_iterations, index;
	double *local_data, *new, *temporary, *send_buffer, *receive_buffer, start_time, stream_time[num_runs+1], exchange_time[num_runs+1], average_time = 0, deviation_time = 0;

	size_t halo_bytes = halo_size*sizeof(double);	// Size of halo region. Will be used later for various memcpys.

	MPI_Comm stream_comm = MPI_COMM_WORLD;		// Communicator for all streaming communications
	MPI_Comm exchange_comm;						// Communicator for all exchange communications
	MPI_Comm_dup(stream_comm,&exchange_comm);

	MPI_Comm_rank(stream_comm, &rank);		// Assigns rank to each processor
	MPI_Comm_size(stream_comm, &size);		// Number of processors involved
	MPI_Request request_send[send_buffer_size], request_receive[receive_buffer_size], send_request_left, send_request_right, receive_request_left, receive_request_right;	// Request arrays are used for streaming, the rest for halo exchange.
	MPI_Status status_send[send_buffer_size], status_receive[receive_buffer_size], send_status_left, send_status_right, receive_status_left, receive_status_right;		// Status arrays are used for streaming, the rest for halo exchange.

	// Divides the array size between the processes if necessary
	if(per_process == 0)
	{
		array_size = array_size/size;
	}

	local_data = malloc((array_size+halo_size)*sizeof(double));		// Primary work buffer. Has enough space for the starting data plus one halo region. Halo streaming stores its work data at the beginning of this array, with the halo region at the end reserved for incoming data. Halo exchange stores its work data in the centre of the array, with space at the beginning and end for halo data.
	new = malloc((array_size+halo_size)*sizeof(double));	// Secondary work buffer to be used in halo exchange.
	send_buffer = malloc((halo_size*message_size*send_buffer_size)*sizeof(double));		// Send buffer for halo streaming
	receive_buffer = malloc((halo_size*message_size*receive_buffer_size)*sizeof(double));	// Receive buffer for halo streaming

	// Parameters used for tracking buffer usage
	int buffer_tracking_index = 0; 		// Index into buffer_usage array
	int array_of_indices[send_buffer_size]; // Indices of outstanding messages
	int buffer_usage[iterations]; 		// Tracks number of outstanding messages
	int max_outstanding = 0;		// Maximum number of outstaqnding messages found
	int *gathered_max_outstanding;
	//int **all_testsome_data;
	int all_testsome_data[iterations][size];
	if (rank==0)
	{
		gathered_max_outstanding = malloc((size)*sizeof(int));
		/*all_testsome_data = malloc((iterations)*sizeof(int*));
		for(i=0;i<iterations;i++)
		{
			all_testsome_data[i] = malloc((size)*sizeof(int));
		}*/
	}


	// Master processor writes parameters to timing file and prints Warnings if there are problems with any of the variables
	if(rank == 0)
	{
		f = fopen(output_filename, "w");
		fprintf(f,"Array size per process: %i\t\tIterations: %i\t\tMessage size: %i\t\tProcesses used: %i\n",array_size,iterations,message_size,size);
		fprintf(f,"Run\t\tStream time\tExchange time\n");

		if(per_process!=0 && per_process!=1)
		{
			printf("WARNING: Non-binary value selected for variable per_process.\n");
		}

		if(send_buffer_size!=receive_buffer_size)
		{
			printf("WARNING: Sizes of Send and Receive Buffers do not match.\n");
		}

		if((array_size%message_size)!=0)
		{
			printf("WARNING: message_size does not divide evenly into array_size per process.\n");
		}

		if((iterations%message_size)!=0)
		{
			printf("WARNING: message_size does not divide evenly into iterations.\n");
		}

		if(message_size>(array_size/2))
		{
			printf("WARNING: message_size is larger than 0.5*array_size. Not enough information to generate a complete message.\n");
		}
	}

	// Loop over number of trials. The +1 means an extra trial is run. This is done so that the time of the first run can be discarded, as it often exhibits cache effects that the other trials do not.
	for(l=0;l<num_runs+1;l++)
	{
		location_send = array_size - halo_size;		// Denotes the point in the array past which the processor no longer has enough information to perform updates. Initialised to array_size - halo_size because there will be halo_size number of points that depend on data from the neighbouring processor.
		location_receive = array_size - halo_size; // Denotes the start of the area that can be calculated when new data arrives. Initialised to the same point as location_send.

		send_iterations = 0;	// Number of iterations worth of data that have been sent
		receive_iterations = 0;	// Number of iterations worth of data that have been received

		// Number of iterations that can be completed without communication. This is either the maximum number of iterations or when the process will run out of data, whichever is smaller.
		int triangle_iterations = (int)fmin((double)array_size/halo_size,(double)iterations);

		// Find neighbours (periodic boundary conditions)
		int neighbour_left = (rank + size - 1) % size;
		int neighbour_right = (rank + 1) % size;

		// Position in global array. This tracks how far the domain of a given processor shifts due to halo streaming.
		int start_point = rank*array_size;

		// Initialise local data
		initialise_data(local_data, array_size);

		// Initialise persistent sends, one for each send tag. These all point to different areas of the Send Buffer. These sends cannot be started as the data is not ready yet.
		for(send_tag=0;send_tag<send_buffer_size;send_tag++)
		{
			MPI_Send_init(&send_buffer[halo_size*message_size*send_tag],halo_size*message_size,MPI_DOUBLE,neighbour_left,send_tag,stream_comm,&request_send[send_tag]);	
		}

		// Initialise persistent receives, one for each receive tag. These all point to different areas of the Receive Buffer. These receives can be started straight away.
		for(receive_tag=0;receive_tag<receive_buffer_size;receive_tag++)
		{
			MPI_Recv_init(&receive_buffer[halo_size*message_size*receive_tag],halo_size*message_size,MPI_DOUBLE,neighbour_right,receive_tag,stream_comm,&request_receive[receive_tag]);
			MPI_Start(&request_receive[receive_tag]);
		}

		// Reset tags so that processing will begin at the correct messages
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

			outstanding_receives(receive_buffer_size,request_receive,buffer_usage,&buffer_tracking_index,array_of_indices,status_receive,&max_outstanding);

			// Wait for request to become available
			MPI_Wait(&request_send[send_tag],&status_send[send_tag]);

			// Completes iterations until a block of data to be sent is complete	
			for(j=0;j<message_size;j++)
			{
				// Copy data from start of work buffer to send buffer 
				memcpy(&send_buffer[halo_size*message_size*send_tag+j*halo_size],local_data,halo_bytes);

				// Update all possible values (in place)
				for(i=0;i<location_send;i++)
				{
					local_data[i] = (local_data[i] + local_data[i+1] + local_data[i+2])/3;
				}

				// Reduce the number of values that can be updated during the next iteration due to less information being available
				location_send -= halo_size;
			}

			// Start new send to transfer the data from the block of iterations that has just been completed
			MPI_Start(&request_send[send_tag]);
			
			// Change tag to next send
			send_tag = (send_tag + 1) % send_buffer_size;

			// Update the number of iterations that have been sent
			send_iterations += message_size;

			// Check if receive has completed
			MPI_Test(&request_receive[receive_tag],&flag,&status_receive[receive_tag]);
		
			if(flag)
			{
				// Contines until all data from the message has been processed
				for(j=0;j<message_size;j++)
				{
					// Copy data from receive buffer to the halo region at the end of the work buffer
					memcpy(&local_data[array_size],&receive_buffer[halo_size*message_size*receive_tag+j*halo_size],halo_bytes);	

					// Update all values that depend on the new data. Stops when it catches up with the data that has been processed for sending
					for(i=0;i+receive_iterations<send_iterations;i++)
					{
						// Cycle over all data in the halo region
						for(k=0;k<halo_size;k++)
						{
							index = array_size + k - halo_size*(i+1);
							local_data[index] = (local_data[index] + local_data[index+1] + local_data[index+2])/3;
						}
					}
					// Update the number of iterations that have been received
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

				// Begin new receive since the data from the last one has all been processed
				MPI_Start(&request_receive[receive_tag]);

				// Change tag to next receive
				receive_tag = (receive_tag + 1) % receive_buffer_size;
			}
		}

		// Extend triangle until its peak reaches the maximum number of iterations
		while(send_iterations<iterations)
		{
			outstanding_receives(receive_buffer_size,request_receive,buffer_usage,&buffer_tracking_index,array_of_indices,status_receive,&max_outstanding);

			// Wait for requests
			MPI_Wait(&request_receive[receive_tag],&status_receive[receive_tag]);
			MPI_Wait(&request_send[send_tag],&status_send[send_tag]);

			// Loop over all received data
			for(j=0;j<message_size;j++)
			{
				// Copy data from receive buffer to the halo region at the end of the work buffer
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

				// Copy data from start of work buffer to send buffer 
				memcpy(&send_buffer[halo_size*message_size*send_tag+j*halo_size],local_data,halo_bytes);
			}

			// Begin new send and receive
			MPI_Start(&request_send[send_tag]);
			MPI_Start(&request_receive[receive_tag]);

			// Move to next set of tags
			receive_tag = (receive_tag + 1) % receive_buffer_size;
			send_tag = (send_tag + 1) % send_buffer_size;

			// Update iterations that have been sent and received
			receive_iterations += message_size;
			send_iterations += message_size;
		}

		// Fill in remaining inverted triangle
		while(receive_iterations<iterations)
		{
			outstanding_receives(receive_buffer_size,request_receive,buffer_usage,&buffer_tracking_index,array_of_indices,status_receive,&max_outstanding);

			// Wait for receive to complete
			MPI_Wait(&request_receive[receive_tag], &status_receive[receive_tag]);

			// Loop over all received data
			for(j=0;j<message_size;j++)
			{
				// Copy data from receive buffer to the halo region at the end of the work buffer
				memcpy(&local_data[array_size],&receive_buffer[halo_size*message_size*receive_tag+j*halo_size],halo_bytes);

				// Update all possible values
				for(i=location_receive;i<array_size;i++)
				{
					local_data[i] = (local_data[i] + local_data[i+1] + local_data[i+2])/3;
				}

				// Increase the number of values to be updated during the next iteration
				location_receive -= halo_size;
			}

			// Begin new receive since the data from the last one has been processed
			MPI_Start(&request_receive[receive_tag]);

			// Change to next receive tag
			receive_tag = (receive_tag + 1) % receive_buffer_size;

			// Update number of iterations received
			receive_iterations += message_size;
		}

		MPI_Barrier(stream_comm);

		// End halo-streaming timing and write time to file
		if (rank == 0)
		{
			stream_time[l] = MPI_Wtime() - start_time;

			if(l==0)
			{
				fprintf(f,"Discarded:\t%f\t",stream_time[l]);
			}
			else
			{
				fprintf(f,"%i\t\t%f\t",l,stream_time[l]);
			}
		}

		// Find new position in global array due to shifting of domains
		start_point = (start_point + iterations) % (size*array_size);

		// Perform I/O operations if it is the last trial
		if(l == num_runs)
		{
			MPI_Gather(&max_outstanding,1,MPI_INT,gathered_max_outstanding,1,MPI_INT,0,stream_comm);
			MPI_Gather(&buffer_usage,iterations,MPI_INT,&all_testsome_data,iterations,MPI_INT,0,stream_comm);
			IO(start_point, array_size, size, local_data, stream_filename);
		}

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
			// Send and receive halo data
			MPI_Isend(&local_data[1],1,MPI_DOUBLE,neighbour_left,j,exchange_comm,&send_request_left);
			MPI_Isend(&local_data[array_size],1,MPI_DOUBLE,neighbour_right,j,exchange_comm,&send_request_right);

			MPI_Irecv(&local_data[0],1,MPI_DOUBLE,neighbour_left,j,exchange_comm,&receive_request_left);
			MPI_Irecv(&local_data[array_size+1],1,MPI_DOUBLE,neighbour_right,j,exchange_comm,&receive_request_right);

			// Update central data
			for(i=2;i<=array_size-1;i++)
			{
				new[i] = (local_data[i-1] + local_data[i] + local_data[i+1])/3;
			}

			// Wait until data has been received
			MPI_Wait(&receive_request_left,&receive_status_left);
			MPI_Wait(&receive_request_right,&receive_status_right);

			// Update edge data
			new[1] = (local_data[0] + local_data[1] + local_data[2])/3;
			new[array_size] = (local_data[array_size-1] + local_data[array_size] + local_data[array_size+1])/3;

			// Wait until data has been sent
			MPI_Wait(&send_request_left,&send_status_left);
			MPI_Wait(&send_request_right,&send_status_right);

			// Swapping pointers to arrays
			temporary = new;
			new = local_data;
			local_data = temporary;
			
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

		// Perform I/O operations if it is the last trial
		if(l == num_runs)
		{

			IO(start_point,array_size,size,&local_data[halo_size/2],halo_filename);
		}
	}

	// Master process computes average and standard deviation of execution times
	if(rank == 0)
	{
		// Compute the average streaming time for all trials except the first one
		for(l=1;l<=num_runs;l++)
		{
			average_time += stream_time[l];
		}

		average_time = average_time/num_runs;

		// Compute the standard deviation of all streaming times except the first one
		for(l=1;l<=num_runs;l++)
		{
			deviation_time += (stream_time[l]-average_time)*(stream_time[l]-average_time);
		}

		deviation_time = sqrt(deviation_time/num_runs);

		// Write the average streaming time to file
		fprintf(f,"Average:\t%f\t",average_time);

		// Reset the average time before beginning calculations for halo exchange
		average_time = 0;

		// Compute the average exchange time for all trials except the first one
		for(l=1;l<=num_runs;l++)
		{
			average_time += exchange_time[l];
		}

		average_time = average_time/num_runs;

		// Write the average exchange time and streaming standard deviation to file
		fprintf(f,"%f\n",average_time);
		fprintf(f,"Deviation:\t%f\t",deviation_time);

		// Reset the standard deviation before beginning calculations on the 
		deviation_time = 0;

		// Compute the standard deviation of all exchange times except the first one
		for(l=1;l<=num_runs;l++)
		{
			deviation_time += (exchange_time[l]-average_time)*(exchange_time[l]-average_time);
		}

		deviation_time = sqrt(deviation_time/num_runs);

		// Write the exchange standard deviation to file
		fprintf(f,"%f\n",deviation_time);
	}

	if(rank==0)
	{
		fprintf(f,"\nMaximum unprocessed messages in buffer for each process:\nProcess\t\tMax Outstanding\n");
		for(i=0;i<size;i++)
		{
			fprintf(f,"%i\t\t%i\n",i,gathered_max_outstanding[i]);
		}

		fprintf(f,"\nAll testsome data\n");

		for(i=0;i<iterations;i++)
		{
			for(j=0;j<size;j++)
			{
				fprintf(f,"%i\t",all_testsome_data[i][j]);
			}
			fprintf(f,"\n");
		}

		/*for(i=0;i<iterations;i++)
		{
    			free(all_testsome_data[i]);
		}

		free(all_testsome_data);*/
		free(gathered_max_outstanding);
	}
	// Free arrays
	free(local_data);
	free(new);
	free(send_buffer);
	free(receive_buffer);

	// Close timing file
	if(rank == 0)
	{
		fclose(f);
	}

	// Finalise MPI
	MPI_Finalize();

	return 0;
}


// Initialise each element of the work buffer to a random integer. The RNG is seeded to make sure the numbers are the same for both the streaming and exchange trials
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
	int leftover = start_point + array_size - array_size*size;	// Used to determine which process contains the data that should go at the start and end of the file
	int i;
	MPI_Status file_status;

	// Only executed by the process that contains the data that should go at the start and end of the file
	if(leftover > 0)
	{
		// Rearrange data before writing to file
		disp = 0;
		double temp[leftover];	// Used to store data temporarily so it can be shifted around in the array

		// Copy the beginning-of-file data to the temporary buffer
		for(i=0;i<leftover;i++)
		{
			temp[i] = local_data[array_size-leftover+i];
		}

		// Shift the end-of-file data to the end of the array
		for(i=array_size-1;i>=leftover;i--)
		{
			local_data[i] = local_data[i - leftover];
		}

		// Copy beginning-of-file data back to the main array
		for(i=0;i<leftover;i++)
		{
			local_data[i] = temp[i];
		}

		// Create split datatype so that data will be written to the correct blocks at the start and end of the file
		int blocklengths[2] = {leftover,array_size-leftover};
		int displacements[2] = {0,(size - 1)*array_size + leftover};
		MPI_Type_indexed(2,blocklengths,displacements,MPI_DOUBLE,&filetype);
	}

	// All other processes creat a simple contiguous datatype representing a block in the middle of the file
	else
	{
		// Displacement into the file
		disp = start_point*sizeof(double);

		MPI_Type_contiguous(array_size,MPI_DOUBLE,&filetype);
	}

	MPI_Type_commit(&filetype);

	MPI_File_open(comm,filename,amode,MPI_INFO_NULL,&fh);

	MPI_File_set_view(fh,disp,MPI_DOUBLE,filetype,"native",MPI_INFO_NULL);

	MPI_File_write_all(fh,local_data,array_size,MPI_DOUBLE,&file_status);

	MPI_File_close(&fh);

	MPI_Type_free(&filetype);
}

void outstanding_receives(int receive_buffer_size, MPI_Request *request_receive, int *buffer_usage, int *buffer_tracking_index, int *array_of_indices, MPI_Status *status_receive, int *max_outstanding)
{
	MPI_Testsome(receive_buffer_size,request_receive,&buffer_usage[*buffer_tracking_index],array_of_indices,status_receive);
	if(buffer_usage[*buffer_tracking_index]>*max_outstanding)
	{
		*max_outstanding = buffer_usage[*buffer_tracking_index];
	}
	*buffer_tracking_index++;
}
