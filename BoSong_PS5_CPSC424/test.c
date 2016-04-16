#include "/home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include "drand.c"
#include <omp.h>


int main(int argc, char* argv){
  int i, j;
  // boundary
  int sum = 0;
  int tmp = 1;
  #pragma omp parallel shared(sum) default(none) private(i) firstprivate(tmp)
  #pragma omp for schedule(static) 
  for(i = 0; i < 100; i++){
    printf("Thread: %d, i = %d, sum = %d\n", omp_get_thread_num(), i, sum);
    sum++;
    
  }
  printf("sum is %d\n", sum);
}
