#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include "headers/io_functions.h"
#include "headers/parameters.h"

// Used for marking the file with the parameters and timestamp
void filename_set(char *filename, char *timeanddate, char *extension, char *num_procs)
{
	sprintf(filename,"%s_AS%i_IT%i_BS%i_MS%i_NP",filename,array_size,iterations,send_buffer_size,message_size);
	strcat(filename,num_procs);
	strcat(filename,timeanddate);
	strcat(filename,extension);
	return;
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

