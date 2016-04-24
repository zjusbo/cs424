#define FP float

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "/home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.h"




void cpu_matrixmult(FP *a,FP *b, FP *c, int n, int p, int m) {

  int index, indexa, indexb;
  FP cvalue;
  for(int col=0;col < m; col++)
    for(int row=0;row < n; row++) {
      indexb = col;
      index = row * m + col;
      cvalue = 0.;
      for (indexa = row*p; indexa < (row*p + p); indexa++, indexb+=m) 
	cvalue += a[indexa]*b[indexb];
      c[index] -= cvalue; //NOTE: This calculates the diff between CPU and GPU computations.
    }
}


int main(int argc, char *argv[]) {

  int i, j; // loop counters

  int gpucount = 0; // Count of available GPUs
  int gpunum = 0; // Device number to use
  int Grid_Dim_x = 1, Grid_Dim_y = 1; //Grid dimension, x and y, square
  int Block_Dim = 1; //Block dimension, x and y, square

  int n, p, m; // matrix dimension
  FP *a,*b,*c;
  FP *dev_a, *dev_b, *dev_c;
  int size_a, size_b, size_c; // number of bytes in arrays

  double elapsed_time_ms; // which is applicable for asynchronous code also
  double etime0, etime1, cptime;
  // --------------------SET PARAMETERS AND DATA -----------------------
  n = atoi(argv[1]);
  p = atoi(argv[2]);
  m = atoi(argv[3]);
   

  size_a = n * p * sizeof(FP); // number of bytes in total in arrays
  size_b = p * m * sizeof(FP);
  size_c = n * m * sizeof(FP);
  a = (FP*) malloc(size_a); // dynamically allocated memory for arrays on host
  b = (FP*) malloc(size_b);
  c = (FP*) malloc(size_c); // results from GPU

  srand(12345);
  // int p = n; //Used here only to illustrate proper initialization for non-square case
  for(i=0;i < n;i++)
    for(j=0;j < p;j++) {
      a[i * p + j] = (FP) rand() / (FP) RAND_MAX;
      //      a[i * p + j] = (FP) i+j; // may be helpful for debugging
    }

  for(i=0;i < p;i++)
    for(j=0;j < m;j++) {
      b[i * m + j] = (FP) rand() / (FP) RAND_MAX;
      //      b[i * n + j] = (FP) i+j; // may be helpful for debugging
    }

  // ------------- COMPUTATION DONE ON GPU ----------------------------

  // cudaEventSynchronize(start); // not needed

timing(&etime0, &cptime);
 cpu_matrixmult(a,b,c, n, p, m); // do calculation on host (NOTE: This computes the diff with GPU result.)
timing(&etime1, &cptime);
elapsed_time_ms = etime1 - etime0;
  printf("Time to calculate results on CPU: %f s.\n", elapsed_time_ms); // exec. time

// ------------------- check device creates correct results -----------------


  free(a);
  free(b);
  free(c);
  return 0;
}
