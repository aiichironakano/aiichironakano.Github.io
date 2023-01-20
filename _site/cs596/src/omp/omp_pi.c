#include <stdio.h>
#include <omp.h>
#define NBIN 100000
#define MAX_THREADS 8

int main() {
	int nthreads,tid;
	double step,sum[MAX_THREADS]={0.0},pi=0.0;

	step = 1.0/NBIN;
#pragma omp parallel private(tid)
	{
		int i;
		double x;
		nthreads = omp_get_num_threads();
		tid = omp_get_thread_num();
		for (i=tid; i<NBIN; i+=nthreads) {
			x = (i+0.5)*step;
			sum[tid] += 4.0/(1.0+x*x);
		}
	}
	for(tid=0; tid<nthreads; tid++)
		pi += sum[tid]*step;
	printf("PI = %f\n",pi);

	return 0;
}
