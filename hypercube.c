#include<stdio.h>
#include<unistd.h>
#include<mpi.h>
#include<math.h>

#define ARR_SIZE 16



int xor(int a, int b){
	return a ^ b;
}

int and(int a,int b){
	return a & b;
}


int main(int argc, char **argv){
double compute_time = 0;
double communication_time = 0;
int c1,p1,b,c,d1,e,i,k,j,p,N1,mask,mask1,dest,dest1,source1,source;  
int world_size;
int rank;
MPI_Status status;
int tag =4;
int master = 0;
int q = 3;
int d = 3;
int chunk_size;
int offset = 0;
int sum = 0;
int value=0;
int partial_sum = 0;
int a[ARR_SIZE];


  MPI_Init(&argc,&argv);
  compute_time = MPI_Wtime();
  MPI_Comm_size(MPI_COMM_WORLD,&world_size);;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  
  chunk_size = ARR_SIZE/world_size;
  offset = chunk_size*rank;

  if(rank == master){
  	for (i=0;i<ARR_SIZE;i++)
  	{
  	a[i]=rand()%10;
  	printf("%d\t",a[i]);
  }}
 mask = (int)pow((double)2,d)-1;
 for(i = d-1;i>=0;i--){
	mask = xor(mask,(int)pow((double)2,i));
	if(and(rank, mask) == 0){
		if(and(rank,(int) pow((double)2,i)) == 0){
			dest = xor(rank,(int) pow((double)2,i));
			communication_time = MPI_Wtime();
			MPI_Send(a,16,MPI_INT,dest,tag,MPI_COMM_WORLD);
		}
		else{
			source = xor(rank,(int)pow((double)2,i));
			MPI_Recv(a,16,MPI_INT,source,tag,MPI_COMM_WORLD,&status);
			communication_time = MPI_Wtime() - communication_time;
		}

	}
	printf("The communcaition time is %lf\n",communication_time);
	}
  mask1 = 0;
    for (j = 0; j<chunk_size; j++)
    {
 	sum = sum + a[j + chunk_size*rank];
    }
    printf("\n SUM IS = %d at Processor %d\n ",sum,rank);
    for (p= 0; p<=d-1; p++)
    {   c1 = pow(2,p);
        if((mask1&rank)==0)
        {
            if((rank&c1)!=0)
            {
                c = pow(2,p);
                dest1 = rank^c;
                MPI_Send(&sum,1,MPI_INT,dest1,0,MPI_COMM_WORLD);
            }
            else
            {
                d1= pow(2,p);
                source1 = rank^d1;
                MPI_Recv(&value,1,MPI_INT,source1,0,MPI_COMM_WORLD,&status);
                sum = sum + value;
            }
         }
	    e= pow(2,p);
            mask1 = mask1^e;
    }

   if(rank==0)
  {
   printf("\nThe total number of processors is %d \n", world_size);
   printf("The sum of the array is %d\n",sum);
   printf("Total Compute tie is %lf\n",MPI_Wtime()-compute_time);
    //time = MPI_Wtime() - time;
   //printf("\nThe total time needed to find the sum: %f\n", time);
  }
MPI_Finalize();
return 0;
}
