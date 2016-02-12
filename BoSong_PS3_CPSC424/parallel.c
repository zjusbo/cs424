#include <stdio.h>
#include <string.h>
#include <stddef.h> 
#include <stdlib.h> 
#include <unistd.h> 
#include "mpi.h"
#include "timing.h"

#define MIN(a,b) (((a)<(b))?(a):(b))

int calBlockLen(int row_col_idx, int block_size);
void block_matmul(double* row, double* col, double* C, int row_idx, int col_idx, int block_size, int N);

int main(int argc, char **argv ) {

  int rank, size, type=99;
  int num_nodes; 
  int worktime;
  double wct0, wct1, wct_comp0, wct_comp1, total_comp_time, cput;
  double wct_comm0, wct_comm1, total_comm_time, total_time;
  if(argc < 2){
    printf("Please provide the size N of matrix.\n");
    return -1;
  } 

  int N, i, j, k;
  double *A, *B, *C;
  int sizeAB, sizeC, iA, iB, iC;
  
  // 1000,2000,4000,8000,12000
  N = atoi(argv[1]); // size of the matrix

  MPI_Status status;

  MPI_Init(&argc,&argv); // Required MPI initialization call

  MPI_Comm_size(MPI_COMM_WORLD,&num_nodes); // Get no. of processes
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Which process am I?
  
  int block_size = N / num_nodes;
  total_comm_time = total_comp_time = 0;
  /* If I am the master (rank 0) ... */
  if (rank == 0) {
    /** Master's work
     * 1. partition A into p rows
     *    assign 2nd ~ pth row to corresponding slave
     * 2. partition B into p columns
     *    assign 2nd ~ pth col to corresponding slave
     * 3. compute result for 1st row and 1st col
     * 4. for i = 1; i < p; i++ 
     *       collect result from node i
     * 5. merge results to build matrix C
     **/
    printf("Master Process: N = %d, num_nodes = %d, rank = %d, block_size = %d.\n", N, num_nodes, rank, block_size);
    sizeAB = N*(N+1)/2; //Only enough space for the nonzero portions of the matrices
    sizeC = N*N; // All of C will be nonzero, in general!

    A = (double *) calloc(sizeAB, sizeof(double)); 
    B = (double *) calloc(sizeAB, sizeof(double)); 
    C = (double *) calloc(sizeC, sizeof(double));
  
    srand(12345); // Use a standard seed value for reproducibility

    // This assumes A is stored by rows, and B is stored by columns. Other storage schemes are permitted
    for (i=0; i<sizeAB; i++) A[i] = ((double) rand()/(double)RAND_MAX);
    for (i=0; i<sizeAB; i++) B[i] = ((double) rand()/(double)RAND_MAX);

    MPI_Barrier(MPI_COMM_WORLD); //wait for everyone to be ready before starting timer
    
    printf("Process %d is behind Barrier now.\n", rank);
    //wct0 = MPI_Wtime(); //set the start time
    timing(&wct_comm0, &cput); //set the start time
    wct0 = wct_comm0;
    /* Send the permanent row to all the workers, which is where the work happens */
    for (i = 1; i < num_nodes; i++) {
      int row_idx = i * block_size;
      iA = row_idx * (row_idx + 1) / 2; // initializes row pointer in A
      // length of flattened block 
      int len = calBlockLen(row_idx, block_size); //gaussian formula
      
      MPI_Send(A + iA, len, MPI_DOUBLE, i, type, MPI_COMM_WORLD); // send permanent row to worker
      printf("Process %d sent row data to process %d.\n", rank, i);
      printf("Process %d: row addr = %p, len = %d, row_idx = %d.\n", rank, A + iA, len, row_idx);
    } 

    /*Send column flow to all workers*/
    for (i = 1; i < num_nodes; i++) {
      int col_idx = i * block_size;
      iB = col_idx * (col_idx + 1) / 2; // initializes col pointer in B
      int len = calBlockLen(col_idx, block_size); //gaussian formula

      MPI_Send(&col_idx, 1, MPI_INT, i, type, MPI_COMM_WORLD);
      printf("Process %d sent to process %d: col_idx = %d.\n", rank, i, col_idx);
      MPI_Send(B + iB, len, MPI_DOUBLE, i, type, MPI_COMM_WORLD);
      printf("Process %d sent to process %d: col addr = %p, len = %d.\n", rank, i, B + iB, len);
    }    
    timing(&wct_comm1, &cput);
    total_comm_time += wct_comm1 - wct_comm0;
    wct_comp0 = wct_comm1;
    //compute results
    double * row = A, * col = B;
    int col_idx = 0;
    int row_idx = 0;
    int row_len = calBlockLen(row_idx, block_size); //gaussian formula
    int col_len = calBlockLen(col_idx, block_size);
    
    block_matmul(row, col, C, row_idx, col_idx, block_size, N);
    printf("Process %d computed first block result\n", rank);

    timing(&wct_comp1, &cput);
    total_comp_time += wct_comp1 - wct_comp0;
    wct_comm0 = wct_comp1;
    //send col to next node, recv col from prev node
    for(i = 1; i < num_nodes; i++){
      
      // send col_idx and col data to next node
      MPI_Send(&col_idx, 1, MPI_INT, rank + 1, type, MPI_COMM_WORLD);
      printf("Process %d sent to process %d: col idx = %d.\n", rank, 1, col_idx);

      MPI_Recv(&col_idx, 1, MPI_INT, num_nodes - 1, type, MPI_COMM_WORLD, &status);
      printf("Process %d recved from process %d: col idx = %d.\n", rank, num_nodes - 1, col_idx);

      MPI_Send(col, col_len, MPI_DOUBLE, rank + 1, type, MPI_COMM_WORLD);
      printf("Process %d sent to process %d: col addr = %p, len = %d.\n", rank, 1, col, col_len);
      
      // recv col_idx and col data from prev node
      col_len = calBlockLen(col_idx, block_size);
      MPI_Recv(col, col_len, MPI_DOUBLE, num_nodes - 1, type, MPI_COMM_WORLD, &status);
      printf("Process %d recved from process %d: col addr = %p, col_len = %d.\n", rank, num_nodes - 1, col, col_len);
      // C is a row block
      timing(&wct_comm1, &cput);
      total_comm_time += wct_comm1 - wct_comm0;
      wct_comp0 = wct_comm1;
      block_matmul(row, col, C, row_idx, col_idx, block_size, N);
      timing(&wct_comp1, &cput);
      total_comp_time += wct_comp1 - wct_comp0;
      wct_comm0 = wct_comp1;
    }
    // collect results from workers
    for(i = 1; i < num_nodes; i++){
      int row_idx = block_size * i;
      MPI_Recv(C + N * row_idx, block_size * N, MPI_DOUBLE, i, type, MPI_COMM_WORLD, &status);
      printf("Process %d recved from process %d: res addr = %p, len = %d.\n", rank, i, C + N * row_idx, block_size * N);
    }
    timing(&wct_comm1, &cput);
    total_comm_time += wct_comm1 - wct_comm0;
    wct1 = wct_comm1;
    total_time = wct1 - wct0;
    printf("Message printed by master: Total elapsed time: %f s. Communication time: %f s. Computation time: %f s\n", total_time, total_comm_time, total_comp_time);

    free(A);
    free(B);
    free(C);    
  }

  /* Otherwise, if I am a worker ... */
  else {
    /**
     * Node i's work
     * 1. recv row from node 0
     * 2. recv col from node 0
     * 3. compute result
     * 4. for j = 1; j < p; j++
     *       send col to node i+1 or 0
     *       recv col from node i-1 or p-1
     *       compute result
     * 5. merge result
     * 6. send result to node 0
     *
     **/
    printf("Slave Process: N = %d, num_nodes = %d, rank = %d, block_size = %d.\n", N, num_nodes, rank, block_size);
    int row_len;
    int col_idx, col_len;
    int row_idx = block_size * rank;
    row_len = calBlockLen(row_idx, block_size);
    double * buf;
    sizeAB = N*(N+1)/2; //Only enough space for the nonzero portions of the matrices
    sizeC = N*N; // All of C will be nonzero, in general!

    A = (double *) calloc(sizeAB, sizeof(double)); 
    B = (double *) calloc(sizeAB, sizeof(double)); 
    C = (double *) calloc(sizeC, sizeof(double));
    buf = (double *) calloc(sizeAB, sizeof(double));

    MPI_Barrier(MPI_COMM_WORLD); //wait for everyone to be ready before starting
    printf("Process %d is behind Barrier now.\n", rank);
    /* Receive permanent row from the master */
    
    printf("Process %d recved from process %d: N = %d.\n", rank, 0, N);
    MPI_Recv(A, row_len, MPI_DOUBLE, 0, type, MPI_COMM_WORLD, &status); // recv permanent row from master
    printf("Process %d recved from process %d: A = %p, row_len = %d.\n", rank, 0, A, row_len);
    
    /* Receive first column from master*/
    MPI_Recv(&col_idx, 1, MPI_INT, 0, type, MPI_COMM_WORLD, &status);
    printf("Process %d recved from process %d: col_idx = %d.\n", rank, 0, col_idx);
    col_len = calBlockLen(col_idx, block_size);
    MPI_Recv(B, col_len, MPI_DOUBLE, 0, type, MPI_COMM_WORLD, &status);
    printf("Process %d recved from process %d: col_addr = %p, col_len = %d.\n", rank, 0, B, col_len);
  
    // calculate result, row_idx pretend to be 0. 
    block_matmul(A, B, C, row_idx, col_idx, block_size, N);
    printf("Process %d: Computer first block completed.\n", rank);
    
    for(i = 1; i < num_nodes; i++){
      int next_rank = rank + 1 == num_nodes? 0: rank + 1;
      int prev_rank = rank - 1;
      int new_col_idx, new_col_len;

      MPI_Recv(&new_col_idx, 1, MPI_INT, prev_rank, type, MPI_COMM_WORLD, &status);
      printf("Process %d recv from process %d: col_idx = %d.\n", rank, prev_rank, col_idx);

      MPI_Send(&col_idx, 1, MPI_INT, next_rank, type, MPI_COMM_WORLD);
      printf("Process %d sent to process %d: col_idx = %d.\n", rank, next_rank, col_idx);

      col_idx = new_col_idx;
      new_col_len = calBlockLen(col_idx, block_size);
      MPI_Recv(buf, new_col_len, MPI_DOUBLE, prev_rank, type, MPI_COMM_WORLD, &status); 
      printf("Process %d recv from process %d: col_len = %d.\n", rank, prev_rank, col_len);
      
      MPI_Send(B, col_len, MPI_DOUBLE, next_rank, type, MPI_COMM_WORLD);
      printf("Process %d sent to process %d: col_len = %d.\n", rank, next_rank, col_len);

      col_len = new_col_len;
      memcpy(B, buf, sizeof(double) * col_len);

      block_matmul(A, B, C, row_idx, col_idx, block_size, N);
      printf("Process %d: block result computed.\n", rank);
    }
    // Send result back to master node
    MPI_Send(C, N * block_size, MPI_DOUBLE, 0, type, MPI_COMM_WORLD);
    printf("Process %d sent to process %d: result sent.\n", rank, 0);
    
    free(A);
    free(B);
    free(C);
    free(buf);
  }

  MPI_Finalize(); // Required MPI termination call
}

int calBlockLen(int row_col_idx, int block_size){
  return (row_col_idx + 1 + row_col_idx + 1 + block_size - 1) * block_size / 2;
}

/**
 * Calculate block multiplication row * col and store the result in C
 * 
 **/
void block_matmul(double* A, double* B, double* C, int row_idx, int col_idx, int block_size, int N){
  int iC, i, iA, iB, j, k;
  for (i = 0; i < block_size; i++) {
    iC = i * N + col_idx;
    iA = i * (row_idx + 1 + (row_idx + i)) / 2; // initializes row pointer in A
    for (j = 0; j < block_size; j++, iC++) {
      iB = j * (col_idx + 1 + (col_idx + j)) / 2; // initializes column pointer in B
      C[iC] = 0.;
      for (k = 0; k <= MIN(i + row_idx, j + col_idx); k++) C[iC] += A[iA+k] * B[iB+k]; // avoids using known-0 entries 
    }
  }
}
