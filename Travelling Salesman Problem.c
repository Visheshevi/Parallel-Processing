/*********************************************
Maram Sulimani
ECE 566
Projject 04
*********************************************/

#include "mpi.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

/**********************************************
 MPI program to solve the Travelling Salesman Problem
 using parallel branch-and-bound tree search algorithm
 *********************************************/
void ringCreation(MPI_Comm *comm1, int *rank, int *world_size);
float *matrixGeneration(int world_size, int rank, int N);
float *ringGeneration(MPI_Comm *comm1, int rank, int world_size, int *N_per_node, int N);

int *localCoordinates;
FILE *file_to_be_read;

int main (int argc, char **argv)
{
   int world_size;
   int rank;
   int i;
   int N;
   int N_per_node;
   float *localMatrix;
   float det[1], result[1];
   double time_for_ring;
   MPI_Status status;
   MPI_Comm comm1;

   file_to_be_read = fopen("dataset.txt", "r");
   char line[250];
   //printf("%s",file_to_be_read);

   fgets(line, 250, file_to_be_read); //discard the first line
   fgets(line, 250, file_to_be_read); //discard 2nd line
   fgets(line, 250, file_to_be_read); //discard 3rd line
   fgets(line, 250, file_to_be_read); //get the line for dimension
   line[strcspn(line, "\r\n")] = '\0'; //strip EPL(s) char at end
   char *token;
   token = strtok(line, " ");
   token = strtok(line, " ");
   token = strtok(line, " ");
   N = atoi(token); //store the number of cities to N
   printf("%d\n",N);
   fgets(line, 250, file_to_be_read); //discard the first line
   fgets(line, 250, file_to_be_read); //discard 2nd line


   MPI_Init(&argc, &argv); //Initialize MPI environment
   MPI_Comm_size(MPI_COMM_WORLD, &world_size); //get the total number of processor

   printf("Size of matrix used is  %d x %d\n", N, N);
   ringCreation(&comm1, &rank, &world_size);
   time_for_ring = MPI_Wtime();
   localMatrix = ringGeneration(&comm1, rank, world_size, &N_per_node,N);


   time_for_ring = MPI_Wtime() - time_for_ring; //total execution time
   printf("\nThe total time needed %f\n", time_for_ring);



   free(localMatrix);
   MPI_Comm_free(&comm1);
   MPI_Finalize();

   return 0;
}

//Ring topology
void ringCreation(MPI_Comm *comm1, int *rank, int *world_size)
{
	MPI_Comm_size(MPI_COMM_WORLD, world_size);

	int dims[1], periods[1];
	dims[0] = *world_size;
	periods[0] = 1;
	localCoordinates = (int*) malloc(sizeof(int) * 1);

	MPI_Cart_create(MPI_COMM_WORLD, 1, dims, periods, 0, comm1);
	MPI_Comm_rank(*comm1, rank);
	MPI_Comm_size(*comm1, world_size);
	MPI_Cart_coords(*comm1, *rank, 1, localCoordinates);
}


//Initialize the matrix with random values {-1, 0 , 1}
float *matrixGeneration(int world_size, int rank, int N)
{
	int i, j;
	float *gen_matrix;
	int low = -1;
	int high = 1;
	double Ts;

	Ts = MPI_Wtime();
	gen_matrix = (float*) malloc(sizeof(float) * (N*N));
	for (i=0; i<N; i++)
	{
		for (j=0; j<N; j++)
		{
			gen_matrix[i*N+j] = (int)(rand() % (high - low + 1) + low);
		}
	}
	Ts = MPI_Wtime() - Ts;
	printf("\nStartup time = %d\n", Ts);
	return gen_matrix;
}

//send the matrix
float *ringGeneration(MPI_Comm *comm1, int rank, int world_size, int *N_per_node,int N)
{
	float *localMatrix;
	int i, j;

	if (rank == 0)
	{
		localMatrix = matrixGeneration(world_size, rank, N);
	}
	else
	{
		localMatrix = (float*) malloc(sizeof(float) * N * N);
	}
	MPI_Bcast(localMatrix, N*N, MPI_FLOAT, 0, *comm1);
	return localMatrix;
}
