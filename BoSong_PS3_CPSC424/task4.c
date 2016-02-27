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
void swap(double** a, double** b);
int cal_block_size(int N, int rank, int num_nodes);
// Use non-blocking function calls

int * _block_size = NULL;
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
  
  int block_size;
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
    //printf("Master Process: N = %d, num_nodes = %d, rank = %d, block_size = %d.\n", N, num_nodes, rank, block_size);
    MPI_Request dummy_request, send_request[20], send1_request[20], recv_request[20];
    MPI_Status status;
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
    //printf("Process %d is behind Barrier now.\n", rank);
    //wct0 = MPI_Wtime(); //set the start time
    timing(&wct_comm0, &cput); //set the start time
    wct0 = wct_comm0;
    /* Send the permanent row to all the workers, which is where the work happens */
    int row_idx = 0;
    for (i = 1; i < num_nodes; i++) {
      block_size = cal_block_size(N, i - 1, num_nodes);
      row_idx += block_size;
      iA = row_idx * (row_idx + 1) / 2; // initializes row pointer in A
      // length of flattened block
      block_size = cal_block_size(N, i, num_nodes); 
      int len = calBlockLen(row_idx, block_size); //gaussian formula
      
      MPI_Isend(A + iA, len, MPI_DOUBLE, i, type, MPI_COMM_WORLD, &send_request[i]); // send permanent row to worker
      //printf("Process %d sent row data to process %d.\n", rank, i);
      //printf("Process %d: row addr = %p, len = %d, row_idx = %d.\n", rank, A + iA, len, row_idx);
    } 
    for(i = 1; i < num_nodes; i++){
      MPI_Wait(&send_request[i], &status);
    }

    int col_idx = 0;
    /*Send column flow to all workers*/
    for (i = 1; i < num_nodes; i++) {
      block_size = cal_block_size(N, i - 1, num_nodes);
      col_idx += block_size;
      iB = col_idx * (col_idx + 1) / 2; // initializes col pointer in B
      block_size = cal_block_size(N, i, num_nodes); 
      int len = calBlockLen(col_idx, block_size); //gaussian formula

      MPI_Isend(&col_idx, 1, MPI_INT, i, type, MPI_COMM_WORLD, &send_request[i]);
      //printf("Process %d sent to process %d: col_idx = %d.\n", rank, i, col_idx);
      MPI_Isend(B + iB, len, MPI_DOUBLE, i, type, MPI_COMM_WORLD, &send1_request[i]);
      //printf("Process %d sent to process %d: col addr = %p, len = %d.\n", rank, i, B + iB, len);
    }    
    for(i = 1; i < num_nodes; i++){
      MPI_Wait(&send_request[i], &status);
    }    
    for(i = 1; i < num_nodes; i++){
      MPI_Wait(&send1_request[i], &status);
    }    

    timing(&wct_comm1, &cput);
    total_comm_time += wct_comm1 - wct_comm0;
    wct_comp0 = wct_comm1;
    //compute results
    double * row = A, * col = B;
    col_idx = 0;
    row_idx = 0;
    int new_col_len, new_col_idx; // used for non blocking
    double * new_col = (double *) calloc(sizeAB, sizeof(double)); 
    double * f_new_col = new_col; // for free use
    block_size = cal_block_size(N, 0, num_nodes);
    int row_len = calBlockLen(row_idx, block_size); //gaussian formula
    int col_len = calBlockLen(col_idx, block_size);
    
    block_matmul(row, col, C, row_idx, col_idx, block_size, N);
    //printf("Process %d computed first block result\n", rank);

    timing(&wct_comp1, &cput);
    total_comp_time += wct_comp1 - wct_comp0;
    wct_comm0 = wct_comp1;
    //send col to next node, recv col from prev node
    for(i = 1; i < num_nodes; i++){
      // send col_idx and col data to next node
      MPI_Isend(&col_idx, 1, MPI_INT, rank + 1, type, MPI_COMM_WORLD, &send_request[0]);
      //printf("Process %d sent to process %d: col idx = %d.\n", rank, 1, col_idx);

      MPI_Irecv(&new_col_idx, 1, MPI_INT, num_nodes - 1, type, MPI_COMM_WORLD, &recv_request[0]);
      //printf("Process %d recved from process %d: col idx = %d.\n", rank, num_nodes - 1, col_idx);

      MPI_Isend(col, col_len, MPI_DOUBLE, rank + 1, type, MPI_COMM_WORLD, &send_request[1]);
     // printf("Process %d sent to process %d: col addr = %p, len = %d.\n", rank, 1, col, col_len);
      
      // recv col_idx and col data from prev node
      
      // block to recv new_col_idx
      MPI_Wait(&recv_request[0], &status); // wait for new_col_idx
      //new_col_len = calBlockLen(new_col_idx, block_size);
      // We can not change this to non-block version
      MPI_Recv(new_col, sizeAB, MPI_DOUBLE, num_nodes - 1, type, MPI_COMM_WORLD, &status);
      //printf("Process %d recved from process %d: col addr = %p, col_len = %d.\n", rank, num_nodes - 1, col, col_len);
      // C is a row block
      MPI_Get_count(&status, MPI_DOUBLE, &new_col_len);
      timing(&wct_comm1, &cput);
      total_comm_time += wct_comm1 - wct_comm0;
      wct_comp0 = wct_comm1;
      block_matmul(row, new_col, C, row_idx, new_col_idx, block_size, N);
      timing(&wct_comp1, &cput);
      total_comp_time += wct_comp1 - wct_comp0;
      wct_comm0 = wct_comp1;
      
      
      MPI_Wait(&send_request[0], &status); // wait to use col idx
      MPI_Wait(&send_request[1], &status); // wait to use col
      col_idx = new_col_idx;
      col_len = new_col_len;
      swap(&col, &new_col);
    }
    // collect results from workers
    row_idx = 0;
    for(i = 1; i < num_nodes; i++){
      block_size = cal_block_size(N, i - 1, num_nodes);  
      row_idx += block_size;
      block_size = cal_block_size(N, i, num_nodes);  
      MPI_Irecv(C + N * row_idx, block_size * N, MPI_DOUBLE, i, type, MPI_COMM_WORLD, &recv_request[i]);
      //printf("Process %d recved from process %d: res addr = %p, len = %d.\n", rank, i, C + N * row_idx, block_size * N);
    }
    for(i = 1; i < num_nodes; i++){
      MPI_Wait(&recv_request[i], &status);
    }    
    timing(&wct_comm1, &cput);
    total_comm_time += wct_comm1 - wct_comm0;
    wct1 = wct_comm1;
    total_time = wct1 - wct0;
    printf("[Process %d] t_total = %fs. t_comm = %fs. t_comp = %fs\n", rank, total_time, total_comm_time, total_comp_time);
    free(A);
    free(B);
    free(C);    
    free(f_new_col);
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
    //printf("Slave Process: N = %d, num_nodes = %d, rank = %d, block_size = %d.\n", N, num_nodes, rank, block_size);
    int row_len;
    int col_idx, col_len;
    int row_idx = 0;
    block_size = cal_block_size(N, rank, num_nodes);
    for(i = 0; i < rank; i++){
      row_idx += cal_block_size(N, i, num_nodes);
    }
    row_len = calBlockLen(row_idx, block_size);
    double * buf;
    sizeAB = N*(N+1)/2; //Only enough space for the nonzero portions of the matrices
    sizeC = N*N; // All of C will be nonzero, in general!
    MPI_Request send_request[20], recv_request[20];
    MPI_Status status;
    A = (double *) calloc(sizeAB, sizeof(double)); 
    B = (double *) calloc(sizeAB, sizeof(double)); 
    C = (double *) calloc(sizeC, sizeof(double));
    buf = (double *) calloc(sizeAB, sizeof(double));

    MPI_Barrier(MPI_COMM_WORLD); //wait for everyone to be ready before starting
    //printf("Process %d is behind Barrier now.\n", rank);
    /* Receive permanent row from the master */
    timing(&wct_comm0, &cput); //set the start time
    
    //printf("Process %d recved from process %d: N = %d.\n", rank, 0, N);
    MPI_Irecv(A, row_len, MPI_DOUBLE, 0, type, MPI_COMM_WORLD, &recv_request[0]); // recv permanent row from master
    //printf("Process %d recved from process %d: A = %p, row_len = %d.\n", rank, 0, A, row_len);
    
    /* Receive first column from master*/
    MPI_Irecv(&col_idx, 1, MPI_INT, 0, type, MPI_COMM_WORLD, &recv_request[1]);
    //printf("Process %d recved from process %d: col_idx = %d.\n", rank, 0, col_idx);
    
    // wait for col_idx
    MPI_Wait(&recv_request[1], &status);
    
    //col_len = calBlockLen(col_idx, block_size);
    
    MPI_Irecv(B, sizeAB, MPI_DOUBLE, 0, type, MPI_COMM_WORLD, &recv_request[2]);
    

    MPI_Wait(&recv_request[0], &status);
    MPI_Wait(&recv_request[2], &status);
    MPI_Get_count(&status, MPI_DOUBLE, &col_len); // length of current col stored in B

    //printf("Process %d recved from process %d: col_addr = %p, col_len = %d.\n", rank, 0, B, col_len);
    timing(&wct_comm1, &cput); //set the start time
    total_comm_time += wct_comm1 - wct_comm0;
    wct_comp0 = wct_comm1;
    // calculate result, row_idx pretend to be 0. 
    block_matmul(A, B, C, row_idx, col_idx, block_size, N);
    //printf("Process %d: Computer first block completed.\n", rank);
    timing(&wct_comp1, &cput); //set the start time
    total_comp_time += wct_comp1 - wct_comp0;
    wct_comm0 = wct_comp1;

    for(i = 1; i < num_nodes; i++){
      int next_rank = rank + 1 == num_nodes? 0: rank + 1;
      int prev_rank = rank - 1;
      int new_col_idx, new_col_len;
      MPI_Irecv(&new_col_idx, 1, MPI_INT, prev_rank, type, MPI_COMM_WORLD, &recv_request[0]);
      //printf("Process %d recv from process %d: col_idx = %d.\n", rank, prev_rank, col_idx);

      MPI_Isend(&col_idx, 1, MPI_INT, next_rank, type, MPI_COMM_WORLD, &send_request[0]);
      //printf("Process %d sent to process %d: col_idx = %d.\n", rank, next_rank, col_idx);

      MPI_Wait(&recv_request[0], &status); // read new_col_idx
      
      //new_col_len = calBlockLen(new_col_idx, block_size);
      MPI_Irecv(buf, sizeAB, MPI_DOUBLE, prev_rank, type, MPI_COMM_WORLD, &recv_request[1]); 
      //printf("Process %d recv from process %d: col_len = %d.\n", rank, prev_rank, col_len);
      
      MPI_Isend(B, col_len, MPI_DOUBLE, next_rank, type, MPI_COMM_WORLD, &send_request[1]);
      //printf("Process %d sent to process %d: col_len = %d.\n", rank, next_rank, col_len);
      
      MPI_Wait(&recv_request[1], &status); // wait for buf
      MPI_Get_count(&status, MPI_DOUBLE, &new_col_len);
      timing(&wct_comm1, &cput);
      total_comm_time += wct_comm1 - wct_comm0;
      wct_comp0 = wct_comm1;
      
      block_matmul(A, buf, C, row_idx, new_col_idx, block_size, N);
      
      timing(&wct_comp1, &cput);
      total_comp_time += wct_comp1 - wct_comp0;
      wct_comm0 = wct_comp1;
      
      MPI_Wait(&send_request[0], &status); // wait for col_idx
      MPI_Wait(&send_request[1], &status); // wait for B
      col_idx = new_col_idx;
      col_len = new_col_len;
      swap(&B, &buf);
    }
    wct_comm0 = wct_comp1;
    // Send result back to master node

    block_size = cal_block_size(N, rank, num_nodes);

    MPI_Send(C, N * block_size, MPI_DOUBLE, 0, type, MPI_COMM_WORLD);
    //printf("Process %d sent to process %d: result sent.\n", rank, 0);
    timing(&wct_comm1, &cput);
    total_comm_time += wct_comm1 - wct_comm0;

    printf("[Process %d] t_comm = %fs, t_comp = %fs\n", rank, total_comm_time, total_comp_time);
    
    free(A);
    free(B);
    free(C);
    free(buf);
  }

  free(_block_size);

  MPI_Finalize(); // Required MPI termination call
}
void swap(double** col, double** new_col){
  double * tmp = *col;
  *col = *new_col;
  *new_col = tmp;
}

// load balance
int cal_block_size(int N, int rank, int num_nodes){
  if(_block_size != NULL){
    // return previously calculated value directly
    return _block_size[rank];
  }
  else{
    // calculate block size for all processes
    int total_len = (1 + N) * N / 2;
    int average_len = total_len / num_nodes;
    _block_size = (int* )malloc(sizeof(int) * num_nodes);
    int node_idx = 0;
    int i;
    int sum = 0;
    int row_idx = 0;
    for(i = 1; i <= N; i++){
      sum += i;
      if(sum > average_len){ // not the last node
        i--;
        _block_size[node_idx++] = i - row_idx;
        row_idx = i;
        sum = 0;
      }
      if(node_idx == num_nodes - 1){
        // last node
        _block_size[node_idx] = N - row_idx;
        break;
      }
    }
    return cal_block_size(N, rank, num_nodes);
  }
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
