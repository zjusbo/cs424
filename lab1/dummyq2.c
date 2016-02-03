#include <stdio.h>

int dummy(double* a){
  printf("this is a dummy function.\n");
  a[0] = 3;
  return 0;
}
