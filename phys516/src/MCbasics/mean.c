/* Monte Carlo integration of PI by sample mean */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
int main() {
  double x, pi, sum = 0.0;
  int try, ntry;
  printf("Input the number of MC trials\n");
  scanf("%d",&ntry);
  srand((unsigned)time((long *)0));
  for (try=0; try<ntry; try++) {
    x = rand()/(double)RAND_MAX;
    sum += 4.0/(1.0 + x*x);
  }
  pi = sum/ntry;
  printf("MC estimate for PI = %f\n", pi);
  return 0;
}
