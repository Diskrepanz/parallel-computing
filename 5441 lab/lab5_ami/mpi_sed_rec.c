#include <mpi.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
	int rank;
	float buffer[5];
	MPI_Status status;

	for(int i =  0; i<5;i++){
		buffer[i] = i;
	}

	MPI_Init( &argc, &argv);
	MPI_Comm_rank( MPI_COMM_WORLD, &rank);// process 0 sending to process 1
	if ( rank == 0 ){//set buffer = something interesting
		MPI_Send( buffer, 5, MPI_FLOAT, 1, 0, MPI_COMM_WORLD);
	}else if ( rank == 1 ){
		MPI_Recv( buffer, 5, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &status);
		printf("Received %d\n", buffer);
	}
	MPI_Finalize();
	return  0;
}