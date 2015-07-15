#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

int main(int argc, char *argv[])
{
	unsigned int array_size = 100000000;
	int iterations = 10;

	double *local_data, *new, *temp, time;
	int i, j, k, num_runs = 10;
	clock_t start, end;

	FILE *f;
	f = fopen("serialdata.dat", "w");

	local_data = malloc(array_size*sizeof(double));
	new = malloc(array_size*sizeof(double));

	for(k=0;k<num_runs;k++)
	{
		srand(227);

		for(i=0;i<array_size;i++)
		{
			local_data[i] = rand()%10000;
		}

		start = clock();
		for(i=0;i<iterations;i++)
		{
			new[0] = (local_data[array_size-1]+local_data[0]+local_data[1])/3;

			for(j=1;j<array_size-1;j++)
			{
				local_data[j] = (local_data[j-1]+local_data[j]+local_data[j+1])/3;
			}

			new[array_size-1] = (local_data[array_size-2]+local_data[array_size-1]+local_data[0])/3;

			temp = new;
			new = local_data;
			local_data = temp;
		}
		end = clock();
		time = ((double)(end - start))/CLOCKS_PER_SEC;
		fprintf(f,"%f\n",time);
	}

	fclose(f);

}
