// Hybrid MPI+CUDA computation of Pi
#include <stdio.h>
#include <mpi.h>
#include <cuda.h>

#define NBIN  10000000  // Number of bins
#define NUM_BLOCK   13  // Number of thread blocks
#define NUM_THREAD 192  // Number of threads per block

// Kernel that executes on the CUDA device
__global__ void cal_pi(float *sum,int nbin,float step,float offset,int nthreads,int nblocks) {
	int i;
	float x;
	int idx = blockIdx.x*blockDim.x+threadIdx.x;  // Sequential thread index across the blocks
	for (i=idx; i<nbin; i+=nthreads*nblocks) {  // Interleaved bin assignment to threads
		x = offset+(i+0.5)*step;
		sum[idx] += 4.0/(1.0+x*x);
	}
}

int main(int argc,char **argv) {
	int myid,nproc,nbin,tid;
	float step,offset,pi=0.0,pig;
	dim3 dimGrid(NUM_BLOCK,1,1);  // Grid dimensions (only use 1D)
	dim3 dimBlock(NUM_THREAD,1,1);  // Block dimensions (only use 1D)
	float *sumHost,*sumDev;  // Pointers to host & device arrays
	int dev_used;

	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);  // My MPI rank
	MPI_Comm_size(MPI_COMM_WORLD,&nproc);  // Number of MPI processes
	nbin = NBIN/nproc;  // Number of bins per MPI process
	step = 1.0/(float)(nbin*nproc);  // Step size with redefined number of bins
	offset = myid*step*nbin;  // Quadrature-point offset

	cudaSetDevice(myid%2);
	size_t size = NUM_BLOCK*NUM_THREAD*sizeof(float);  //Array memory size
	sumHost = (float *)malloc(size);  //  Allocate array on host
	cudaMalloc((void **) &sumDev,size);  // Allocate array on device
	cudaMemset(sumDev,0,size);  // Reset array in device to 0
	// Calculate on device (call CUDA kernel)
	cal_pi <<<dimGrid,dimBlock>>> (sumDev,nbin,step,offset,NUM_THREAD,NUM_BLOCK);
	// Retrieve result from device and store it in host array
	cudaMemcpy(sumHost,sumDev,size,cudaMemcpyDeviceToHost);
	// Reduction over CUDA threads
	for(tid=0; tid<NUM_THREAD*NUM_BLOCK; tid++)
		pi += sumHost[tid];
	pi *= step;
	// CUDA cleanup
	free(sumHost);
	cudaFree(sumDev);
	cudaGetDevice(&dev_used);
	printf("myid = %d: device used = %d; partial pi = %f\n",myid,dev_used,pi);
	// Reduction over MPI processes
	MPI_Allreduce(&pi,&pig,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
	if (myid==0) printf("PI = %f\n",pig);

	MPI_Finalize();
	return 0;
}
