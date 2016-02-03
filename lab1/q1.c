/**
 * Author: Bo Song
 * Title: CSPS 424/524
 * Date: 02/02/2016
 **/
#include <stdio.h>
#include <stdlib.h>

int SLICES;
double delta_x, Pi, x, sum;
void init(void* arg){
  SLICES = 1000000;
  delta_x = 1./SLICES;
  sum = 0.;
}

void run(void* arg){
  for(int i = 0; i < SLICES; i++){
    x = (i + 0.5) * delta_x;
    sum += (1.0 / (1.0 + x * x));
  }
  Pi = 4.0 * sum * delta_x;
}
