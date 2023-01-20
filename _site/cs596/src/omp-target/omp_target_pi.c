#include <omp.h>
#include <stdio.h>
#define NBIN 1000000000
#define NTRD 1024

int main() {
	double step,sum=0.0,pi;
	step = 1.0/(double)NBIN;

	#pragma omp target map(step,sum)
	{
		long long i;
		double x;
		# pragma omp parallel for private (i,x) reduction(+:sum) num_threads(NTRD)
		for (i=0; i<NBIN; i++) {
			x = (i+0.5)*step;
			sum += 4.0/(1.0+x*x);
		}
	}
	pi = sum*step;
	printf("PI = %f\n",pi);
	return 0;
}
