/**
 * Author: Bo Song
 * Title: CSPS 424/524
 * Date: 02/02/2016
 * Description: estimate division delay
 **/
#include <stdio.h>
#include <stdlib.h>

int SLICES;
volatile double sum;
void init(void* arg){
  SLICES = 1000000;
}

void run(void* arg){
  for(int i = 0; i < SLICES; i++){
    sum = 1.0 / i;
  }
  
}
