#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <omp.h>
static uint64_t seed;
#define MAX_NUM 2812500
#pragma omp threadprivate(seed)


void advance(int);
void dsrand(unsigned s)
{
        int i, n;

	seed = s-1;
        i = omp_get_thread_num();
        n = omp_get_num_threads();
        advance(i * (MAX_NUM / n));
        printf("Seed = %lu. RAND_MAX = %d.\n",seed,RAND_MAX);
        
}

double drand(void)
{
  double re;
  {
	seed = 6364136223846793005ULL*seed + 1;
        re = ((double)(seed>>33)/(double)RAND_MAX);
  }
  return re;
}

void advance(int n){
  int i;
  printf("Thread: %d, advance %d\n", omp_get_thread_num(), n);
  for(i = 0; i < n; i++){
    seed = 6364136223846793005ULL*seed + 1;
  }
}
