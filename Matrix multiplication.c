#include<stdio.h>
#include<unistd.h>
#include<mpi.h>

#define k 4
#define size 4


//Method to find the cofactors foro the determinants
void getCofactor(long long int mat[size][size],long long int temp[size][size], int p, int q, int n){
    int i = 0, j = 0;
    int row,col;
    for (row = 0; row < n; row++){
        for (col = 0; col < n; col++){
            if (row != p && col != q){
                temp[i][j++] = mat[row][col];
 		if (j == n - 1){
                    j = 0;
                    i++;
                }
            }
        }
    }
}

//Method to find the determinant of the final matrix
long long int determinantOfMatrix(long long int mat[size][size], int n){
  long long int D = 0;
  int f=0;
  if (n == 1)
        return mat[0][0];
  long long int temp[size][size];
  int sign = 1;
  for(f = 0;f<n;f++){
	getCofactor(mat,temp,0,f,n);
	D += sign*mat[0][f] * determinantOfMatrix(temp,n-1);
	sign = -sign;
  }
  return D;
}
int main(int argc, char **argv){
  int world_size = 0;
  int rank = 0;
  long long int mat_a[size*size];
  long long int mat_b[size*size];
  long long int mat_final[size][1];
  long long int temp_1[size];
  long long int temp_2[size];
  int i=0,j=0;
  long long int recv[size*size];
  int tag = 4;
  long long int temp_mat[size*size];
  int num = 0;
  int big_loop=0;
  long long int fin_mat[size][size];
  long long int determinant=0;
  double initial_time_mul=0;
  double final_time_mul=0;
  double initial_computation_time = 0;
  double final_computation_time = 0;
  double initial_send_time;
  double final_send_time;
  int high = 1;
  int low = -1;
  MPI_Status status;


  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&world_size); //Find world size
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);  //Find the rank of the processor
  if(rank == 0){
  	initial_time_mul = MPI_Wtime();
  	initial_computation_time = MPI_Wtime();
  }
  if(rank == 0){
	//Initialize the matrix in processor at rank 0
	long long int scope_a[size*size];
	long long int scope_b[size*size];

	for(i=0;i<size*size;i++){
		scope_a[i] = (long long int)(rand() % (high-low+1)+low);
		scope_b[i] = scope_a[i];
	}
	for(i=0;i<size*size;i++){
		printf("%lld\t",scope_a[i]);
	}
	printf("\n");
	for(i = 0;i<size*size;i++){
		mat_a[i] = scope_a[i];
		mat_b[i] = scope_b[i];
	}
  }
  for(big_loop = 0;big_loop<k-1;big_loop++){
  	if(rank ==0){
		//Send the matrix to other processors
		for(i = 1;i<world_size;i++){
        		initial_send_time = MPI_Wtime();
			MPI_Send(&mat_a[0],size*size,MPI_LONG_LONG_INT,i,tag,MPI_COMM_WORLD);
			MPI_Send(&mat_b[0],size*size,MPI_LONG_LONG_INT,i,tag,MPI_COMM_WORLD);
			final_send_time = MPI_Wtime() - initial_send_time;
		}
 	 }
	//Receive the matrix from the processor at rank 0
  	if(rank!=0){
		MPI_Recv(&mat_a[0],size*size,MPI_LONG_LONG_INT,0,tag,MPI_COMM_WORLD,&status);
		MPI_Recv(&mat_b[0],size*size,MPI_LONG_LONG_INT,0,tag,MPI_COMM_WORLD,&status);
  	}

	//Extract the required data from the matrix
  	for(i = rank;i<size*size;i=i+size){
		temp_1[num] = mat_a[i];
  		num++;
  	}
  	num = 0;
  	for(i = rank*size;i<(rank*size)+size;i++){
		temp_2[num] = mat_b[i];
		num++;
  	}
  	num = 0;

 	//Calculate the intermediate matrix for eeach processor
  	for(i = 0;i<size;i++){
		for(j = 0;j<size;j++){
			temp_mat[num] = temp_1[i]*temp_2[j];
			num++;
		}
  	}
        num = 0;
	
	//Sending hte intermediate matrix from all processors to processor at rank 0 and then adding them
	if(rank!=0){
		MPI_Send(&temp_mat[0],size*size,MPI_LONG_LONG_INT,0,tag,MPI_COMM_WORLD);
	}
	if(rank == 0){
		for(i = 1;i<world_size;i++){
			MPI_Recv(&recv[0],size*size,MPI_LONG_LONG_INT,i,tag,MPI_COMM_WORLD,&status);
			for(j=0;j<size*size;j++){
				temp_mat[j] = temp_mat[j] + recv[j];
			}
		}
	}

  	if(rank ==0){
  		for(i = 0;i<size*size;i++){
			mat_a[i] = temp_mat[i];
  		}
		num = 0;
		for(i=0;i<size;i++){
			for(j=0;j<size;j++){
				fin_mat[i][j] = mat_a[num];
				num++;
			}
		}
		num = 0;
  	}
  }
  if(rank ==0)
  	final_time_mul = MPI_Wtime() - initial_time_mul;

  if(rank ==0){
	for(i = 0;i<size*size;i++)
		printf("%lld\t",mat_a[i]);
	printf("\n");
  }
  if(rank ==0){
	double ts = MPI_Wtime();
	determinant = determinantOfMatrix(fin_mat,size);
	double fts = MPI_Wtime() - ts;
	final_computation_time = MPI_Wtime() - initial_computation_time;
	printf("The determinant of matrix is %lld\n",determinant);
	printf("Time for determinant is %lf\n",fts);
	printf("Time to multply is %lf\n",final_time_mul);
	printf("Total Compuation Time %lf\n",final_computation_time);
  	printf("Send time is %lf\n",final_send_time);
  }
  MPI_Finalize();

/*
  for(int i = 0;i<3;i++){
	for(int j = 0;j<3;j++){
		printf("%d\t",mat_a[i][j]);
		}
	printf("\n");
  }
*/
/*
 if(rank ==0){
	printf("%d\t%d\t%d\n",mat_b[0][0],mat_b[1][0],mat_b[2][0]);
  }
*/

  return 0;
}
