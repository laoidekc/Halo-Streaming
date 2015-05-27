#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
	MPI_Init(NULL, NULL);


	// Parameter declarations
	int array_size = 1000000;
	int iterations = 1000;
	int send_buffer_size = 10;
	int halo_size = 2;

	int rank, size, i, j, location = array_size-halo_size, tag, leftover;
	double *local_data, *new, *temporary, *send_buffer, start_time, time;

	local_data = malloc((array_size+halo_size)*sizeof(double));
	new = malloc((array_size+halo_size)*sizeof(double));
	send_buffer = malloc((halo_size*send_buffer_size)*sizeof(double));

	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &size);
	MPI_Request request[send_buffer_size], request_left, request_right;
	MPI_Status status_send[send_buffer_size], status_receive[send_buffer_size], status_left, status_right, file_status;
	int amode = MPI_MODE_CREATE | MPI_MODE_WRONLY;
	MPI_File fh;
	MPI_Datatype filetype;
	MPI_Offset disp;

	// Finding neighbours
	int neighbour_left = (rank + size - 1) % size;
	int neighbour_right = (rank + 1) % size;

	// Position in global array
	int start_point = rank*array_size;

	// Initialising local data
	for(i=0;i<array_size;i++)
	{
		local_data[i] = (double)i;
	}

	// Initialising persistent sends
	for(tag=0;tag<send_buffer_size;tag++)
	{
		MPI_Send_init(&send_buffer[halo_size*tag],halo_size,MPI_DOUBLE,neighbour_left,tag,comm,&request[tag]);	
	}

	tag = 0;

	MPI_Barrier(comm);

	// Beginning halo-streaming timing
	if(rank == 0)
	{
		start_time = MPI_Wtime();
	}

	for(j=0;j<iterations;j++)
	{
		// Waiting for request to become available
		MPI_Wait(&request[tag],&status_send[tag]);

		// Copying data to send buffer 
		for(i=0;i<halo_size;i++)
		{
			send_buffer[halo_size*tag+i] = local_data[i];
		}

		// Sending data
		MPI_Start(&request[tag]);

		// Changing tag
		tag = (tag + 1) % send_buffer_size;

		// Updating all possible values (in place)
		for(i=0;i<location;i++)
		{
			local_data[i] = (local_data[i] + local_data[i+1] + local_data[i+2])/3;
		}

		// Reducing the number of values that can be updated during the next iteration
		location -= halo_size;
	}
	
	// Waiting for all sends to complete (this may be unnecessary)
	MPI_Waitall(send_buffer_size,request,status_send);

	// Finding new position in global array
	start_point = (start_point + iterations) % (size*array_size);

	// Finding the point where new data will be received
	location = array_size - halo_size;
	tag = 0;

	for(j=0;j<iterations;j++)
	{
		// Receiving data from neighbour
		MPI_Recv(&local_data[array_size],halo_size,MPI_DOUBLE,neighbour_right,tag,comm,&status_receive[tag]);

		// Changing tag
		tag = (tag + 1) % send_buffer_size;

		// Updating in place based on new data
		for(i=location;i<array_size;i++)
		{
			local_data[i] = (local_data[i] + local_data[i+1] + local_data[i+2])/3;
		}

		// Changing the location where received data is stored
		location -= halo_size;
	}

	MPI_Barrier(comm);

	// Ending halo-streaming timing
	if (rank == 0)
	{
		time = MPI_Wtime() - start_time;
		printf("Halo streaming time taken = %lf seconds\n",time);
	}

	// Opening file
	MPI_File_open(comm,"out.bin",amode,MPI_INFO_NULL,&fh);

	// Determing which process must write to both ends of the file
	leftover = start_point + array_size - array_size*size;

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

	// Setting file view, writing to file, and closing file
	MPI_File_set_view(fh,disp,MPI_DOUBLE,filetype,"native",MPI_INFO_NULL);

	MPI_File_write_all(fh,local_data,array_size,MPI_DOUBLE,&file_status);

	MPI_File_close(&fh);

	MPI_Type_free(&filetype);


	// Re-initialising data for halo-exchange
	for(i=1;i<=array_size;i++)
	{
		local_data[i] = (double)i-1;
	}

	MPI_Barrier(comm);

	// Beginning timing for halo-exchange
	if(rank == 0)
	{
		start_time = MPI_Wtime();
	}

	for(j=0;j<iterations;j++)
	{
		// Sending and receiving halos
		MPI_Issend(&local_data[1],1,MPI_DOUBLE,neighbour_left,tag,comm,&request_left);
		MPI_Issend(&local_data[array_size],1,MPI_DOUBLE,neighbour_right,tag,comm,&request_right);

		MPI_Recv(&local_data[array_size+1],1,MPI_DOUBLE,neighbour_right,tag,comm,&status_right);
		MPI_Recv(&local_data[0],1,MPI_DOUBLE,neighbour_left,tag,comm,&status_left);

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

	MPI_Barrier(comm);

	// Ending halo-exchange timing
	if (rank == 0)
	{
		time = MPI_Wtime() - start_time;
		printf("Halo exchange time taken = %lf seconds\n",time);
	}

	// Opening file
	MPI_File_open(comm,"halo_out.bin",amode,MPI_INFO_NULL,&fh);
	disp = rank*array_size*sizeof(double);

	// Creating contiguous datatype
	MPI_Type_contiguous(array_size,MPI_DOUBLE,&filetype);
	MPI_Type_commit(&filetype);

	// Setting file view, writing to file, and closing file
	MPI_File_set_view(fh,disp,MPI_DOUBLE,filetype,"native",MPI_INFO_NULL);
	MPI_File_write_all(fh,&local_data[1],array_size,MPI_DOUBLE,&file_status);
	MPI_File_close(&fh);

	// Freeing arrays
	free(local_data);
	free(new);
	free(send_buffer);

	MPI_Finalize();

	return;
}
