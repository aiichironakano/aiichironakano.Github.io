#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main() {
  int Max_step; /* Maximum number of random-walk steps */
  int N_walker; /* Number of walkers */
  int step,walker,k;
  int x; /* Drunkard's position */
  int hist[1001];

  /* Input parameters */
  printf("Input the number of maximum random-walk steps\n");
  scanf("%d",&Max_step);
  printf("Input the number of walkers\n");
  scanf("%d",&N_walker);

  for (k=0; k<=1000; k++)
    hist[k] = 0.0;

  srand((unsigned)time((long *)0)); /* Initialize the rondom-number sequence */

  for (walker=1; walker<=N_walker; walker++) {
    x = 0;
    for (step=1; step<=Max_step; step++) {
      if (rand() < RAND_MAX/2)
        ++x;
      else
        --x;
    } /* Endfor step */
    k = x + 500;
    ++hist[k];
  } /* Endfor walker */

  for (k=0; k<1001; k++)
    printf("%d %d\n",k-500,hist[k]);

  return 0;
}
