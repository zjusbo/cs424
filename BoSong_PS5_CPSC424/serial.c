#include "/home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include "drand.c"
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
double stride = 0.001;

int main(int argc, char* argv){
  Compl z, c, tmp;
  double i, j;
  dsrand(12345);
  int N1, N0;
  double etime0, etime1, cptime;
  double A; 
  int n;
  // boundary
  lowerLeft.real = -2.0;
  lowerLeft.imag = 0;
  upperRight.real = 0.5;
  upperRight.imag = 1.125;
 
  N1 = N0 = 0;
  timing(&etime0, &cptime);
  for(i = lowerLeft.real; i < upperRight.real; i += stride){
    for(j = lowerLeft.imag; j < upperRight.imag; j += stride){
      c.real = drand() * stride + i;
      c.imag = drand() * stride + j;
      z = c;
      for(n = 0; n < MAX_ITERATE; n++){
        multiply(&z, &c, &tmp);
        z = tmp;
        if(lengthsq(&z) > threshold * threshold){
            break;
        }        
      }
      if(n == MAX_ITERATE){
        N1++;
      }else{
        N0++;
      }
    }
  }
  timing(&etime1, &cptime);
  A = 2.0 * N1 / (N0 + N1) * (upperRight.real - lowerLeft.real) * (upperRight.imag - lowerLeft.imag);
  printf("area is %f, time elapsed is %f\n", A, etime1 - etime0);
}
