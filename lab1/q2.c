#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>



extern int dummyq2(double*);
extern void timing(double*, double*);

double * a, * b, * c, * d;
int n; 
time_t t;
double randd();
double mypow(double base, int exp){
  double re = 1;
  for(int i = 0; i < exp; i++){
    re *= base;
  }
  return re; 
}
void init(void * arg){
  int p;
  long long ll;
  double * pd = (double*) &ll;
  p = atoi(((char **)arg)[1]);
  n = (int)mypow(2.1, p); 

  srand((unsigned)time(&t));

  a = (double *)calloc(n, sizeof(double));
  b = (double *)calloc(n, sizeof(double));
  c = (double *)calloc(n, sizeof(double));
  d = (double *)calloc(n, sizeof(double));

  for(int i = 0; i < n; i++){
    a[i] = randd();
    b[i] = randd();
    c[i] = randd();
    d[i] = randd();
  }
} 

double randd(){
  long long ll;
  ll = ((long long)rand() << 32) + rand();
  return *(double *)&ll;
}



/**
 * Author: Bo Song
 * Title: CPSC424 Lab1
 * Date: 02/02/2016
 **/

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
      for(int i = 0; i < n; i++){
        a[i] = b[i] + c[i] * d[i];
      }
      if (t < 1000) dummyq2(a); // foo the compiler
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
