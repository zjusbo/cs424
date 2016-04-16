#include "/home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <omp.h>
#include <stdint.h>


extern void dsrand(unsigned int);
extern double drand();

typedef struct complextype {
  double real, imag;
} Compl;

void multiply(const Compl* z, const Compl* c, Compl* result){
  result->real = z->real * z->real - z->imag * z->imag + c->real;
  result->imag = 2.0 * z->real * z->imag + c->imag;  
}

double lengthsq(const Compl* z){
  return z->real * z->real + z->imag * z->imag;
}
Compl lowerLeft, upperRight;
int MAX_ITERATE = 20000;
double threshold = 2.0;
double stride = 0.001; // it should be 0.001. Here I use 0.1 for quick test

int main(int argc, char* argv){
  uint64_t seed;
  Compl z, c, tmp;
  int i, j;
  int N1, N0;
  double etime0, etime1, cptime;
  double A; 
  double r;
  int n;
  // boundary
  lowerLeft.real = -2.0;
  lowerLeft.imag = 0;
  upperRight.real = 0.5;
  upperRight.imag = 1.125;
  omp_set_dynamic(0);
  #pragma omp parallel
  dsrand(12345);

  N1 = N0 = 0;
  timing(&etime0, &cptime);

  #pragma omp parallel firstprivate(j, n, c, z, tmp, stride) \
                      shared(i, N0, N1, lowerLeft, upperRight, threshold, MAX_ITERATE)
  #pragma omp for schedule(static, 10) collapse(2)
  for(i = 0; i < (int)((upperRight.real - lowerLeft.real) / stride); i++){
    for(j = 0; j < (int)((upperRight.imag - lowerLeft.imag) / stride); j++){
      if(i == 0 && j == 0 && omp_get_thread_num() == 0) printf("Threads: %d\n", omp_get_num_threads());
      c.real = lowerLeft.real + (drand() + i) * stride;
      c.imag = lowerLeft.imag + (drand() + j) * stride;
      z = c;
      for(n = 0; n < MAX_ITERATE; n++){
        multiply(&z, &c, &tmp);
        z = tmp;
        if(lengthsq(&z) > threshold * threshold){
            break;
        }        
      }
      if(n == MAX_ITERATE){
        #pragma omp critical(N1)
        N1++;
      }else{
        #pragma omp critical(N0)
        N0++;
      }
    }
  }
  timing(&etime1, &cptime);
  A = 2.0 * N1 / (N0 + N1) * (upperRight.real - lowerLeft.real) * (upperRight.imag - lowerLeft.imag);
  printf("area is %f, time elapsed is %f, total cell: %d\n", A, etime1 - etime0, N1 + N0);
}
