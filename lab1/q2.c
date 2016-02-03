#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

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

void run(void* arg){
  for(int i = 0; i < n; i++){
    a[i] = b[i] + c[i] * d[i];
  }
}
double randd(){
  long long ll;
  ll = ((long long)rand() << 32) + rand();
  return *(double *)&ll;
}
