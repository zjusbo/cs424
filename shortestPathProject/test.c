#include <omp.h>
#include <stdio.h>

int main(int argc, char ** argv){
  #pragma omp parallel
  printf("%d\n", omp_get_num_threads());

  return 0;
}
