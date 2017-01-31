#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "parameters.h"

void initialise_data(double*);
void IO(int, int, double*, char*); 
void filename_set(char*, char*, char*, char*);

int main(int argc, char *argv[])
{
	MPI_Init(NULL, NULL);

	int rank, size;

	MPI_Comm stream_comm = MPI_COMM_WORLD;		// Communicator for all streaming communications
	MPI_Comm exchange_comm;						// Communicator for all exchange communications
	MPI_Comm_dup(stream_comm,&exchange_comm);	// Communicators are identical and only exist to separate communications

	MPI_Comm_rank(stream_comm, &rank);		// Assigns rank to each processor
	MPI_Comm_size(stream_comm, &size);		// Number of processors involved

	// parameters.c allows the choice of whether to run the streaming code, exchange code, or both. This makes sure at least one is selected and exits if not
	if(rank==0 && !run_streaming && !run_exchange)
	{
		printf("Error: Must run either exchange code or streaming code\n");
		MPI_Finalize();
		return 0;
	}

	// Output files
	char stream_filename[256] = "streaming_out";		// Binary output of the final halo streaming data
	char halo_filename[256] = "exchange_out";	// Binary output of the final halo exchange data
	char output_filename[256] = "data";			// Will contain the timing results

	// Gets current time in string. This is used for timestamping the output files
	time_t current_time;
	struct tm * time_info;
	char time_string[21];
	time(&current_time);
	time_info = localtime(&current_time);
	strftime(time_string,sizeof(time_string),"_%F_%H-%M-%S",time_info);

	// Adds timestamp and extension to filenames
	filename_set(stream_filename,time_string,".bin",argv[1]);
	filename_set(halo_filename,time_string,".bin",argv[1]);
	filename_set(output_filename,time_string,".txt",argv[1]);
	FILE *f;

	int i, j, k, l, location_send, location_receive, send_tag, receive_tag, leftover, flag, send_iterations, receive_iterations, index;
	double *local_data, *new, *temporary, *send_buffer, *receive_buffer, start_time, stream_time[num_runs+1], exchange_time[num_runs+1], average_time = 0, deviation_time = 0;

	size_t halo_bytes = halo_size*sizeof(double);	// Size of halo region. Will be used later for various memcpys.

	MPI_Request request_send[send_buffer_size], request_receive[receive_buffer_size], send_request_left, send_request_right, receive_request_left, receive_request_right;	// Request arrays are used for streaming, the rest for halo exchange.
	MPI_Status status_send[send_buffer_size], status_receive[receive_buffer_size], send_status_left, send_status_right, receive_status_left, receive_status_right;		// Status arrays are used for streaming, the rest for halo exchange.

	// Divides the array size between the processes if necessary. Set in parameters.c
	if(per_process == 0)
	{
		array_size = array_size/size;
	}

	local_data = malloc((array_size+halo_size)*sizeof(double));		// Primary work buffer. Has enough space for the starting data plus one halo region. Halo streaming stores its work data at the beginning of this array, with the halo region at the end reserved for incoming data. Halo exchange stores its work data in the centre of the array, with space at the beginning and end for halo data.
	new = malloc((array_size+halo_size)*sizeof(double));	// Secondary work buffer to be used in halo exchange.
	send_buffer = malloc((halo_size*message_size*send_buffer_size)*sizeof(double));		// Send buffer for halo streaming
	receive_buffer = malloc((halo_size*message_size*receive_buffer_size)*sizeof(double));	// Receive buffer for halo streaming

	// Find neighbours (periodic boundary conditions)
	int neighbour_left = (rank + size - 1) % size;
	int neighbour_right = (rank + 1) % size;

	// Used to track where the local data is located in the global dataset
	int start_point;


	// Master processor writes parameters to timing file and prints warnings if there are problems with any of them
	if(rank == 0)
	{
		f = fopen(output_filename, "w");
		fprintf(f,"#Array size per process: %i\t\tIterations: %i\t\tMessage size: %i\t\tProcesses used: %i\n#Run\t\t",array_size,iterations,message_size,size);
		if(run_streaming)
		{
			fprintf(f,"Stream time\t");
		}		

		if(run_exchange)
		{
			fprintf(f,"Exchange time");
		}

		fprintf(f,"\n");

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
		// Write run number to timing file
		if(rank==0)
		{
			if(l==0)
			{
				fprintf(f,"#Discarded:\t");
			}
			else
			{
				fprintf(f,"%i\t\t",l);
			}
		}

		// Halo streaming portion
		if(run_streaming)
		{
			location_send = array_size - halo_size;		// Denotes the point in the array past which the processor no longer has enough information to perform updates. Initialised to array_size - halo_size because there will be halo_size number of points that depend on data from the neighbouring processor.
			location_receive = array_size - halo_size; // Denotes the start of the area that can be calculated when new data arrives. Initialised to the same point as location_send.

			send_iterations = 0;	// Number of iterations worth of data that have been sent
			receive_iterations = 0;	// Number of iterations worth of data that have been received

			// Number of iterations that can be completed without communication. This is either the maximum number of iterations or when the process will run out of data, whichever is smaller.
			int triangle_iterations = (int)fmin((double)array_size/halo_size,(double)iterations);

			// Initialises position of local data in global dataset
			start_point = rank*array_size;

			// Initialise local data
			initialise_data(local_data);

			// Initialise persistent sends, one for each send tag. These all point to different areas of the Send Buffer. These sends cannot be started as the data is not ready yet.
			for(send_tag=0;send_tag<send_buffer_size;send_tag++)
			{
				MPI_Send_init(&send_buffer[halo_size*message_size*send_tag], // There are send_buffer_size number of Sends initialised. These are separated by halo_size*message_size in memory, which is the size of each message
								halo_size*message_size,                      // Size of message
								MPI_DOUBLE,                                  // MPI Datatype
								neighbour_left,                              // Receiver
								send_tag,                                    // Message tag
								stream_comm,                                 // Communicator
								&request_send[send_tag]);                    // Request
			}

			// Initialise persistent receives, one for each receive tag. These all point to different areas of the Receive Buffer. Arguments are the same as those of the Sends above
			for(receive_tag=0;receive_tag<receive_buffer_size;receive_tag++)
			{
				MPI_Recv_init(&receive_buffer[halo_size*message_size*receive_tag],
								halo_size*message_size,
								MPI_DOUBLE,
								neighbour_right,
								receive_tag,
								stream_comm,
								&request_receive[receive_tag]);
				MPI_Start(&request_receive[receive_tag]); // Can start Receives immediately
			}

			// Reset tags so that processing will begin at the correct messages
			send_tag = 0;
			receive_tag = 0;

			// Synchronise so timing can begin
			MPI_Barrier(stream_comm);

			// Begin halo-streaming timing
			if(rank == 0)
			{
				start_time = MPI_Wtime();
			}

			// Compute initial triangle
			while(send_iterations<triangle_iterations)
			{
				// Wait for send request to become available. This should return instantly unless ALL Sends are outstanding
				MPI_Wait(&request_send[send_tag],
								&status_send[send_tag]);

				// Completes iterations until a block of data to be sent is complete
				for(j=0;j<message_size;j++)
				{
					// Copy data from start of work buffer to send buffer 
					memcpy(&send_buffer[halo_size*message_size*send_tag+j*halo_size], // halo_size*message_size*send_tag takes you to the start of the current message. halo_size*the loop index takes you to ythe slot for the data from the current iteration
									local_data, // Start of work buffer
									halo_bytes); // Message size in bytes

					// Update every piece of data for which the necessary data is available (in place)
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
				MPI_Test(&request_receive[receive_tag],
								&flag,
								&status_receive[receive_tag]);
				
				if(flag)
				{
					// Contines until all data from the message has been processed
					for(j=0;j<message_size;j++)
					{
						// Copy data from receive buffer to the halo region at the end of the work buffer. Indexing works the same as in the Send buffer case.
						memcpy(&local_data[array_size],
										&receive_buffer[halo_size*message_size*receive_tag+j*halo_size],
										halo_bytes);	

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

					// Increases amount of data that can be calculated in the new iteration, due to extra information being available
					location_send += halo_size*message_size;

					// Increase size of triangle that can be computed, due to new information having arrived
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
				// Wait for requests
				MPI_Wait(&request_receive[receive_tag],
								&status_receive[receive_tag]);
				MPI_Wait(&request_send[send_tag],
								&status_send[send_tag]);

				// Loop over all received data
				for(j=0;j<message_size;j++)
				{
					// Copy data from receive buffer to the halo region at the end of the work buffer.
					memcpy(&local_data[array_size],
									&receive_buffer[halo_size*message_size*receive_tag+j*halo_size],
									halo_bytes);

					// Calculate all values that depend on the new data
					for(i=array_size-halo_size;i>=0;i-=halo_size)
					{
						// Cycle over all data in the halo region
						for(k=0;k<halo_size;k++)
						{
							index = i+k;
							local_data[index] = (local_data[index] + local_data[index+1] + local_data[index+2])/3;
						}
					}

					// Copy data from start of work buffer to send buffer.
					memcpy(&send_buffer[halo_size*message_size*send_tag+j*halo_size],
									local_data,
									halo_bytes);
				}

				// Begin new send and receive
				MPI_Start(&request_send[send_tag]);
				MPI_Start(&request_receive[receive_tag]);

				// Move to next set of tags
				receive_tag = (receive_tag + 1) % receive_buffer_size;
				send_tag = (send_tag + 1) % send_buffer_size;

				// Update number of iterations that have been sent and received
				receive_iterations += message_size;
				send_iterations += message_size;
			}

			// Fill in remaining inverted triangle
			while(receive_iterations<iterations)
			{
				// Wait for receive to complete. (Hopefully this should never have to wait)
				MPI_Wait(&request_receive[receive_tag],
								&status_receive[receive_tag]);

				// Loop over all received data
				for(j=0;j<message_size;j++)
				{
					// Copy data from receive buffer to the halo region at the end of the work buffer.
					memcpy(&local_data[array_size],
									&receive_buffer[halo_size*message_size*receive_tag+j*halo_size],
									halo_bytes);

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

			// Ensure all processes have finished
			MPI_Barrier(stream_comm);

			// End halo-streaming timing and write time to file
			if (rank == 0)
			{
				stream_time[l] = MPI_Wtime() - start_time;
				fprintf(f,"%f\t",stream_time[l]);
			}

			// Find new position in global array due to shifting of domains
			start_point = (start_point + iterations) % (size*array_size);

			// Perform I/O operations if it is the last trial
			if(l == num_runs)
			{
				IO(start_point, size, local_data, stream_filename);
			}

			// Cancel all outstanding communications
			for(receive_tag=0;receive_tag<receive_buffer_size;receive_tag++)
			{
				MPI_Cancel(&request_receive[receive_tag]);
			}
		}

		if(run_exchange)
		{
			// Re-initialise data for halo-exchange
			initialise_data(&local_data[halo_size/2]);

			// Synchronise processes for timing purposes
			MPI_Barrier(exchange_comm);

			// Begin timing for halo-exchange
			if(rank == 0)
			{
				start_time = MPI_Wtime();
			}

			for(j=0;j<iterations;j++)
			{
				// Send and receive halo data
				MPI_Isend(&local_data[1],1,MPI_DOUBLE,neighbour_left,0,exchange_comm,&send_request_left);
				MPI_Isend(&local_data[array_size],1,MPI_DOUBLE,neighbour_right,0,exchange_comm,&send_request_right);

				MPI_Irecv(&local_data[0],1,MPI_DOUBLE,neighbour_left,0,exchange_comm,&receive_request_left);
				MPI_Irecv(&local_data[array_size+1],1,MPI_DOUBLE,neighbour_right,0,exchange_comm,&receive_request_right);

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

			// Ensure all processes have finished
			MPI_Barrier(exchange_comm);

			// End halo-exchange timing and write time to file
			if(rank == 0)
			{
				exchange_time[l] = MPI_Wtime() - start_time;
				fprintf(f,"%f",exchange_time[l]);
			}

			// Find position in global array. This is independent of the number of iterations for halo exchange
			start_point = array_size*rank;

			// Perform I/O operations if it is the last trial
			if(l == num_runs)
			{
				IO(start_point,size,&local_data[halo_size/2],halo_filename);
			}
		}

		// Starts new line in timing file for next run's data
		if(rank==0)
		{
			fprintf(f,"\n");
		}
	}

	// Master process computes average and standard deviation of execution times
	if(rank == 0)
	{
		fprintf(f,"#Average\t");
		if(run_streaming)
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
			fprintf(f,"%f\t",average_time);
		}
	
		if(run_exchange)
		{
			// Reset the average time before beginning calculations for halo exchange
			average_time = 0;

			// Compute the average exchange time for all trials except the first one
			for(l=1;l<=num_runs;l++)
			{
				average_time += exchange_time[l];
			}

			average_time = average_time/num_runs;

			// Write the average exchange time and streaming standard deviation to file
			fprintf(f,"%f",average_time);
		}

		fprintf(f,"\n#Deviation:\t");
		
		if(run_streaming)
		{
			fprintf(f,"%f\t",deviation_time);
		}

		if(run_exchange)
		{
			// Reset the standard deviation before beginning calculations on the 
			deviation_time = 0;

			// Compute the standard deviation of all exchange times except the first one
			for(l=1;l<=num_runs;l++)
			{
				deviation_time += (exchange_time[l]-average_time)*(exchange_time[l]-average_time);
			}

			deviation_time = sqrt(deviation_time/num_runs);

			// Write the exchange standard deviation to file
			fprintf(f,"%f",deviation_time);
		}

		fprintf(f,"\n");
	}

	// Free arrays
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
void initialise_data(double *local_data)
{
	srand(227);
	int i;

	for(i=0;i<array_size;i++)
	{
		local_data[i] = rand()%10000;
	}
}

void IO(int start_point, int size, double *local_data, char *filename)
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

	// All other processes create a simple contiguous datatype representing a block in the middle of the file
	else
	{
		// Displacement into the file
		disp = start_point*sizeof(double);

		MPI_Type_contiguous(array_size,MPI_DOUBLE,&filetype);
	}

	// MPI-IO stuff that does the actual writing to the file
	MPI_Type_commit(&filetype);

	MPI_File_open(comm,filename,amode,MPI_INFO_NULL,&fh);

	MPI_File_set_view(fh,disp,MPI_DOUBLE,filetype,"native",MPI_INFO_NULL);

	MPI_File_write_all(fh,local_data,array_size,MPI_DOUBLE,&file_status);

	MPI_File_close(&fh);

	MPI_Type_free(&filetype);
	return;
}

// Used for marking the file with the parameters and timestamp
void filename_set(char *filename, char *timeanddate, char *extension, char *num_procs)
{
	sprintf(filename,"%s_AS%i_IT%i_BS%i_MS%i_NP",filename,array_size,iterations,send_buffer_size,message_size);
	strcat(filename,num_procs);
	strcat(filename,timeanddate);
	strcat(filename,extension);
	return;
}
