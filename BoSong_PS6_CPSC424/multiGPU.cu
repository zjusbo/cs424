#include <stdio.h>
#include <cuda.h>
#include <stdlib.h>
#include <math.h>


#define FP float
#define TW 32
#define N 4
#define DEVCOUNT 2

__global__ void gpu_matrixmult(FP *a,FP *b, FP *c, int n, int p, int m) {

  __shared__ FP atile[TW][TW], 
             btile[N][TW][TW]; 
  int tx = threadIdx.x; int ty = threadIdx.y; 
  FP cvalue[N];
  int col[N];
  for(int i = 0; i < N; i++){
    col[i] = tx + blockDim.x * (blockIdx.x * N + i);
    cvalue[i] = 0;
  }
  int row = ty + blockDim.y * blockIdx.y;
  int num_tile = N; //(m - blockDim.x * blockIdx.x * N) / TW  < N? (m - blockDim.x * blockIdx.x * N) / TW: N ; 
  for(int k = 0; k < p / TW; k++){
    atile[ty][tx] = a[row * p + k * TW + tx];
    for(int i = 0; i < num_tile; i++)
      btile[i][ty][tx] = b[(k * TW + ty) * m + col[i]];
         
      __syncthreads();
       for(int i = 0; i < num_tile; i++)
        for(int l = 0; l < TW; l++) cvalue[i] += atile[ty][l] * btile[i][l][tx];

      __syncthreads();
  }
    for(int i = 0; i < num_tile; i++) 
      c[row*m + col[i]] = cvalue[i];
}

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
  FP *a,*b,*c_part[DEVCOUNT], *c;
  FP *dev_a[DEVCOUNT], *dev_b[DEVCOUNT], *dev_c[DEVCOUNT];
  int size_a, size_b, size_c; // number of bytes in arrays

  cudaEvent_t start, stop; // using cuda events to measure time
  float elapsed_time_ms; // which is applicable for asynchronous code also
  cudaError_t errorcode;

  // --------------------SET PARAMETERS AND DATA -----------------------

  errorcode = cudaGetDeviceCount(&gpucount);
  if (errorcode == cudaErrorNoDevice) {
    printf("No GPUs are visible\n");
    exit(-1);
  }
  else {
     printf("Device count = %d\n",gpucount);
  }

  if ((argc<5) || (argc>6)) {
    printf("Usage: matmul <matrix dim n> <matrix dim p> <matrix dim m> <block dim> [<dev num>]\n");
    exit (-1);
  }

  n = atoi(argv[1]);
  p = atoi(argv[2]);
  m = atoi(argv[3]);
   
  Block_Dim = atoi(argv[4]); // Square block
  if (Block_Dim*Block_Dim > 1024) {
    printf("Error, too many threads in block\n");
//    exit (-1);
  }
  Grid_Dim_x = (m / Block_Dim  + N - 1) / N; 
  Grid_Dim_y = (n / DEVCOUNT ) / Block_Dim; // divide gride by 2 vertically

  if (argc==6) {
    gpunum = atoi(argv[5]); // Device number
    if ((gpunum > 2) || (gpunum < 0)) {
      printf("Error, Device number must be 0, 1, or 2\n");
      exit (-1);
    }
  }
  printf("Using device %d\n",gpunum);
  
  printf("Matrix Dimension = (%d, %d, %d)\n",n, p, m);
  printf("Block_Dim = %d, Grid_Dim_x = %d, Grid_Dim_y = %d\n",Block_Dim, Grid_Dim_x, Grid_Dim_y);

  dim3 Grid(Grid_Dim_x, Grid_Dim_y); //Grid structure
  dim3 Block(Block_Dim, Block_Dim); //Block structure

  size_a = n * p * sizeof(FP); // number of bytes in total in arrays
  size_b = p * m * sizeof(FP);
  size_c = n * m * sizeof(FP);
  a = (FP*) malloc(size_a); // dynamically allocated memory for arrays on host
  b = (FP*) malloc(size_b);
  c = (FP*) malloc(size_c);
  for(int i = 0; i < DEVCOUNT; i++){
    c_part[i] = (FP*) malloc(size_c / DEVCOUNT);  
  }
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
  // allocate space
  for(int i = 0; i < DEVCOUNT; i++){
    cudaSetDevice(i);
    cudaMalloc((void**)&dev_a[i], size_a / DEVCOUNT); // allocate memory on device
    cudaMalloc((void**)&dev_b[i], size_b);
    cudaMalloc((void**)&dev_c[i], size_c / DEVCOUNT);
   
    cudaMemcpyAsync(dev_a[i], a + i * n * p / DEVCOUNT, size_a / DEVCOUNT ,cudaMemcpyHostToDevice);
    cudaMemcpyAsync(dev_b[i], b , size_b ,cudaMemcpyHostToDevice);

  }

  // run kernel
  cudaSetDevice(0);
  cudaEventCreate(&start); // instrument code to measure start time
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0); 
  for(int i = 0; i < DEVCOUNT; i++){
    cudaSetDevice(i);
    gpu_matrixmult<<<Grid,Block>>>(dev_a[i],dev_b[i],dev_c[i],n / DEVCOUNT, p, m); // compute the left part of matrix
  } 

  // wait for completion
  cudaSetDevice(0);
  cudaEventRecord(stop, 0); // instrument code to measure end time
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elapsed_time_ms, start, stop);

  // copy result back to host
  for(int i = 0; i < DEVCOUNT; i++){
    cudaSetDevice(i);
    cudaMemcpy(c_part[i],dev_c[i], size_c / DEVCOUNT,cudaMemcpyDeviceToHost);
  }

  // merge results
  for(int i = 0; i < DEVCOUNT; i++){
    memcpy(c + m * n * i / DEVCOUNT, c_part[i], size_c / DEVCOUNT);
  }
  printf("Time to calculate results on GPU: %f ms.\n", elapsed_time_ms); // exec. time
  cudaSetDevice(0);
  // ------------- COMPUTATION DONE ON HOST CPU ----------------------------
  // DEBUGGING USE ONLY (AND FOR LIMITED NUMBERS OF TIMING RUNS)

  cudaEventRecord(start, 0); // use same timing
  // cudaEventSynchronize(start); // not needed

  cpu_matrixmult(a,b,c, n, p, m); // do calculation on host (NOTE: This computes the diff with GPU result.)

  cudaEventRecord(stop, 0); // instrument code to measue end time
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elapsed_time_ms, start, stop );

  printf("Time to calculate results on CPU: %f ms.\n", elapsed_time_ms); // exec. time

// ------------------- check device creates correct results -----------------

  double error, suma, sumb, sumc, ai, bi, ci;
  suma = 0.; sumb = 0; sumc = 0;
  for(i=0;i < n*m;i++) {
    ai = (double) a[i];
    bi = (double) b[i];
    ci = (double) c[i];
    suma += ai*ai;
    sumb += bi*bi;
    sumc += ci*ci;
  }
  suma = sqrt(suma);
  sumb = sqrt(sumb);
  sumc = sqrt(sumc);
  error =  sumc/(n*suma*sumb);
  printf("Scaled error between GPU and CPU: %e\n", error);

// -------------- clean up ---------------------------------------

  free(a);
  free(b);
  free(c);
  for(int i = 0; i < DEVCOUNT; i++){
    free(c_part[i]);
    cudaFree(dev_a[i]);
    cudaFree(dev_b[i]);
    cudaFree(dev_c[i]);
  }
    cudaEventDestroy(start);
    cudaEventDestroy(stop);


  return 0;
}
