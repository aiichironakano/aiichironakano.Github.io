#include <omp.h>
#include <stdio.h>
#define NBIN 1000000
#define NTRD 96

int main() {
	float step,sum=0.0,pi;
	step = 1.0/(float)NBIN;

	#pragma omp target map(step,sum)
	{
		# pragma omp parallel for reduction(+:sum) num_threads(NTRD)
		for (long long i=0; i<NBIN; i++) {
			float x = (i+0.5)*step;
			sum += 4.0/(1.0+x*x);
		}
	}

	pi = sum*step;
	printf("PI = %f\n",pi);
	return 0;
}
