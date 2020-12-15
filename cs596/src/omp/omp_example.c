#include <stdio.h>
#include <omp.h>

int main () {
	int nthreads,tid;

	nthreads = omp_get_num_threads();
	printf("Sequential section: # of threads = %d\n",nthreads);

	/* Fork multi-threads with own copies of variable */
	#pragma omp parallel private(tid) 
	{
		/* Obtain & print thread id */
		tid = omp_get_thread_num();
		printf("Parallel section: Hello world from thread %d\n",tid);
		/* Only master thread does this */
		if (tid == 0) {
			nthreads = omp_get_num_threads();
			printf("Parallel section: # of threads = %d\n",nthreads);
		}
	} /* All created threads terminate */

	return 0;
}
