#include <stdio.h>
#include <mpi.h>
#include <omp.h>
#define NBIN 100000
#define MAX_THREADS 8

int main(int argc,char **argv) {
	int nbin,myid,nproc,nthreads,tid;
	double step,sum[MAX_THREADS]={0.0},pi=0.0,pig;

	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nproc);
	nbin = NBIN/nproc;

	step = 1.0/(nbin*nproc);

	omp_set_num_threads(2);

#pragma omp parallel private(tid)
	{
		int i;
		double x;
		nthreads = omp_get_num_threads();
		tid = omp_get_thread_num();
		for (i=nbin*myid+tid; i<nbin*(myid+1); i+=nthreads) {
			x = (i+0.5)*step;
			sum[tid] += 4.0/(1.0+x*x);
		}
		printf("rank tid sum = %d %d %e\n",myid,tid,sum[tid]);
	}
	for (tid=0; tid<nthreads; tid++)
		pi += sum[tid]*step;

	MPI_Allreduce(&pi,&pig,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	if (myid==0) printf("PI = %f\n",pig);

	MPI_Finalize();
	return 0;
}
