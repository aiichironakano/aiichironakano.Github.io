/* Hypercube mergesort using MPI **********************************************/
#include <stdio.h>
#include <math.h>
#include "mpi.h"

#define N 1024 /* Maximum list size */
#define MAX 99 /* Maximum value of a list element */

int nprocs,dim,myid; /* Cube size, dimension, & my rank */

/* Sequential mergesort (either ascending or descending) */
void mergesort(int list[],int left,int right,int descending)
{
  int i,j,k,t,middle,temp[N];

  if (left < right) {
    middle = (left + right)/2;
    mergesort(list, left, middle, descending);
    mergesort(list, middle+1, right, descending);

    k = i = left; j = middle+1;
    if (descending)
      while (i<=middle && j<=right)
        temp[k++] = list[i]>list[j] ? list[i++] : list[j++];
    else
      while (i<=middle && j<=right)
        temp[k++] = list[i]<list[j] ? list[i++] : list[j++];
    t = i>middle ? j : i;
    while (k <= right) temp[k++] = list[t++];
    for (k=left; k<=right; k++) list[k] = temp[k];
  } 
}

/* Parallel mergesort */
void parallel_mergesort(int myid,int list[],int n)
{
  int listsize, l, m, bitl = 1, bitm, partner, i;
  MPI_Status status;

  listsize = n/nprocs;
  mergesort(list,0,listsize-1,myid & bitl);

  for (l=1; l<=dim; l++) {
    bitl = bitl << 1;
    for (bitm=1, m=0; m<l-1; m++) bitm *= 2;
    for (m=l-1; m>=0; m--) {
      partner = myid ^ bitm;
      MPI_Send(list,listsize,MPI_INT,partner,l*dim+m,MPI_COMM_WORLD);
      MPI_Recv(&list[listsize],listsize,MPI_INT,partner,l*dim+m,
               MPI_COMM_WORLD,&status);
      mergesort(list,0,2*listsize-1,myid & bitl);
      if (myid & bitm)
        for (i=0; i<listsize; i++) list[i] = list[i+listsize];
      bitm = bitm >> 1;
    }
  } 
}

int main(int argc, char *argv[])
{
  int list[N],n=16,i;

  MPI_Init(&argc,&argv); /* Initialize the MPI environment */

  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  dim = log(nprocs+1e-10)/log(2.0);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

  srand((unsigned) myid+1);
  for (i=0; i<n/nprocs; i++) list[i] = rand()%MAX;

  printf("Before: Rank %2d :",myid);
  for (i=0; i<n/nprocs; i++) printf("%3d ",list[i]);
  printf("\n");

  parallel_mergesort(myid,list,n);

  printf("After:  Rank %2d :",myid);
  for (i=0; i<n/nprocs; i++) printf("%3d ",list[i]);
  printf("\n");

  MPI_Finalize();

  return 0;
}
