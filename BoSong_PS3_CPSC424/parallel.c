#include <stdio.h>
#include <string.h>
#include <stddef.h> 
#include <stdlib.h> 
#include <unistd.h> 
#include "mpi.h"
#include "timing.h"

main(int argc, char **argv ) {

  /*
    This is the Hello World program for CPSC424/524.

    Author: Andrew Sherman, Yale University

    Date: 9/15/2014

    Credits: This program is based on a program provided by Barry Wilkinson (UNCC), which 
             had a similar communication pattern, but did not include any simulated work.
  */

  char message[100];
  int i,rank, size, type=99; 
  int worktime, sparm, rwork(int,int);
  double wct0, wct1, total_time, cput;

  MPI_Status status;

  MPI_Init(&argc,&argv); // Required MPI initialization call

  MPI_Comm_size(MPI_COMM_WORLD,&size); // Get no. of processes
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Which process am I?
  /* If I am the master (rank 0) ... */
  char msg_buf[size + 1][100];
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
    sparm = rwork(0,0); //initialize the workers' work times
    
    /* Create the message using sprintf */
    sprintf(message, "Hello, from process %d.",rank);

    MPI_Barrier(MPI_COMM_WORLD); //wait for everyone to be ready before starting timer
    //wct0 = MPI_Wtime(); //set the start time
    timing(&wct0, &cput); //set the start time

    /* Send the message to all the workers, which is where the work happens */
    for (i=1; i<size; i++) {
      MPI_Send(message, strlen(message)+1, MPI_CHAR, i, type, MPI_COMM_WORLD);
      MPI_Send(&sparm, 1, MPI_INT, i, type, MPI_COMM_WORLD);
    } 

    //wct1 = MPI_Wtime(); // Get total elapsed time
    for(i = 1; i < size; i++){
      MPI_Recv(message, 100, MPI_CHAR, MPI_ANY_SOURCE, type, MPI_COMM_WORLD, &status);
      strcpy(msg_buf[status.MPI_SOURCE], message);
      sleep(3);
    }
    for(i = 1; i < size; i++){
      printf("Message from process %d: %s\n", i, msg_buf[i]);
    }
    timing(&wct1, &cput); //get the end time
    total_time = wct1 - wct0;
    printf("Message printed by master: Total elapsed time is %f seconds.\n",total_time);

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
    MPI_Barrier(MPI_COMM_WORLD); //wait for everyone to be ready before starting
  /* Receive messages from the master */
    MPI_Recv(message, 100, MPI_CHAR, 0, type, MPI_COMM_WORLD, &status);
    MPI_Recv(&sparm, 1, MPI_INT, 0, type, MPI_COMM_WORLD, &status);

    worktime = rwork(rank,sparm); // Simulate some work
    sprintf(message, "Hello master, from process %d after working %d seconds.",
	   rank, worktime);
    MPI_Send(message, strlen(message) + 1, MPI_CHAR, 0, type, MPI_COMM_WORLD);
  }

  MPI_Finalize(); // Required MPI termination call
}
