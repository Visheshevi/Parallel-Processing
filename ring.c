#include <mpi.h>
#include <stdio.h>

int main(int argc, char **argv) {
  
  int left, right;
  int N_per_process;
   int i, j;
   int start, end;
   double time;    
   int sum = 0;
   int ownSum=0;
   MPI_Status status;
   
  
  int world_size;
  int rank;
  char hostname[256];
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  int master = 0;
  int v=0;
  int tag = 0;
  int N;
  int factor=0;
  int t=0;
  const int size=1000; 
  int array[size];
  
  MPI_Init(&argc,&argv);

 
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

 
  MPI_Get_processor_name(processor_name, &name_len);

   gethostname(hostname,255);
   
   printf("\nThe total number of processors: %d\n", world_size);
   
   N=size/world_size;
   factor=N*rank;
     for(i = 0;i<size;i++){
	array[i] = rand() % 10;
  }	
   for (i=0;i<N;i++)
   	{
   		ownSum=ownSum+array[i+factor];
	   }

   if(rank==0)
   {
   	   double t1=MPI_Wtime();
	   MPI_Send(&sum,1,MPI_INT,(rank+1),tag,MPI_COMM_WORLD);
	   MPI_Recv(&sum,1,MPI_INT,world_size-1,0,MPI_COMM_WORLD,&status);
	   double t2=MPI_Wtime();
	   sum=sum+ownSum;
	   printf("the sum of the array: %d",sum);
	   double tt=t2-t1;
	   printf("the time is: %lf", tt);
   }
   if((rank>0)&&(rank<world_size-1))
   {
   	MPI_Recv(&sum,1,MPI_INT,rank-1,0,MPI_COMM_WORLD,&status);
   	sum=sum+ownSum;
   	MPI_Send(&sum,1,MPI_INT,(rank+1),tag,MPI_COMM_WORLD);
   }
   if(rank==world_size-1)
   {
   	MPI_Recv(&sum,1,MPI_INT,rank-1,0,MPI_COMM_WORLD,&status);
    sum=sum+ownSum;
   	MPI_Send(&sum,1,MPI_INT,0,tag,MPI_COMM_WORLD);
   	}
   
   time = MPI_Wtime() - time; 
   printf("\nThe total time is needed to find the sum: %f\n", time);
   

   
   
  MPI_Finalize();

  return 0;
}
