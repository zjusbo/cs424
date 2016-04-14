#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "/home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.h"
#include <mpi.h>
/*
  This is the serial program for CPSC424/524 Assignment #4.
  
  Author: Andrew Sherman, Yale University
  Modified by: Bo Song, Yale University
  Date: 2/23/2013

*/
// Globals so that we don't need to pass them to functions
double dt, dt2;
int N, n, K;
double *mass, *x, *y, *z, *vx, *vy, *vz, *fx, *fy, *fz, *buf;
double *bmass, *bx, *by, *bz; // boundary bodies
// Main function
int main(int argc, char **argv) {
  int i, ts, thisbody, otherbody;
  int bn;
  double vavgx, vavgy, vavgz, ax, ay, az, deltaf[3];
  // Serial Timing Functions: See the timer.h include file in this directory
  double etime, etime0, etime1, cptime;

  int send_counts[8], recv_counts[8], send_offsets[8], recv_offsets[8];  
  int rank, size, type=99;
  MPI_Init(&argc, &argv);  
  MPI_Comm_size(MPI_COMM_WORLD, &size); // Get no. of processes
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Which process am I?


  MPI_Alltoall(send_counts, 1, MPI_INT, recv_counts, 1, MPI_INT, MPI_COMM_WORLD); 
  printf("Program finish\n");


  MPI_Finalize();
}
