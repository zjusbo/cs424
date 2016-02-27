#include "stdio.h"
#include "mpi.h"
#include "stdlib.h"
#include "stddef.h"

typedef struct car_s{
	int length;
	double weight;
	int color[3];
}car;
int main(int argc, char ** argv){
	const int tag = 13;
	int size, rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	const int nitems = 3;
	// how many elements in each block
	int blocklength[] = {1, 2, 3};
	MPI_Datatype types[] = {MPI_INT, MPI_DOUBLE, MPI_INT};
	MPI_Datatype mpi_car_type;
	MPI_Aint offsets[nitems];

	offsets[0] = offsetof(car, length);
	offsets[1] = offsetof(car, weight);
	offsets[2] = offsetof(car, color);

	MPI_Type_create_struct(nitems, blocklength, offsets, types, &mpi_car_type);
	MPI_Type_commit(&mpi_car_type);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank == 0){
		car send;
		send.length = 10;
		send.weight = 23.4;
		send.color[0] = 1;
                send.color[1] = 2;
                send.color[2] = 3;
		const int dest = 1;
                MPI_Request req;
		MPI_Isend(&send, 1, mpi_car_type, dest, tag, MPI_COMM_WORLD, &req);
		printf("Master finish\n");
	}
	else{
		MPI_Status status;
		car recv;
		const int src = 0;
		MPI_Recv(&recv, 1, mpi_car_type, src, tag, MPI_COMM_WORLD, &status);
		printf("Process %d Recieved: length = %d, weight = %f, color[0] = %d, color[1] = %d, color[2] = %d",
			rank, recv.length, recv.weight, recv.color[0], recv.color[1], recv.color[2]);
	}
	MPI_Type_free(&mpi_car_type);
	MPI_Finalize();

}
