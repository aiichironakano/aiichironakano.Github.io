/* Monte Carlo integration of PI by hit and miss */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
int main() {
  double x, y, pi;
  int inner,hit = 0, try, ntry;
  printf("Input the number of MC trials\n");
  scanf("%d",&ntry);
  srand((unsigned)time((long *)0));
  for (try=0; try<ntry; try++) {
    x = rand()/(double)RAND_MAX;
    y = rand()/(double)RAND_MAX;
    if (x*x+y*y < 1.0) ++hit;
  }
  pi = 4.0*hit/ntry;
  printf("Ntry Hit = %d %d\n",ntry,hit);
  printf("MC estimate for PI = %f\n", pi);
  return 0;
}