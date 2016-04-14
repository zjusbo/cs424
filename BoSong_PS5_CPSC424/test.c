#include <stdio.h>

typedef struct s_A{
 int a;
 int b;

} A;

int main(int argc, char ** argv){
  A a;
  a.a = 1;
  a.b = 2;
  A b;
  b = a;
  b.a = 2;

  printf("%d %d\n", a.a, a.b);
  printf("%d %d\n", b.a, b.b);

  return 0;
}
