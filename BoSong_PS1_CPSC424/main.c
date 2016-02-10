/**
 * Author: Bo Song
 * Title: CPSC424 Lab1
 * Date: 02/02/2016
 **/
#include <stdio.h>
#include <time.h>

extern int dummy();
extern void timing(double*, double*);
extern void run(void *); // kernal benchmark function
extern void init(void *);  // init benchmark environment 

int main(int argc, char ** argv){
  int repeat = 1;
  double runtime = 0.;
  double wcs, wce, ct;
  time_t t;
  time(&t);
  init((void *)argv);
  for(; runtime <.1; repeat *= 2){
    timing(&wcs, &ct);
    for(int r = 0; r < repeat; ++r){
      //kern code
      run((void*)argv);
      if (t < 1000) dummy(); // foo the compiler
    }
    timing(&wce, &ct);
    runtime = wce - wcs;
  }
  repeat /= 2;
  printf("Repeat for %d times. Wallclock time is %f.", repeat, runtime);
  if(argc > 1){
    printf("Argv is %s.", argv[1]);
  }
  printf("\n");
  return 0;
}
