#include "mpi.h"
#include <stdio.h>
#define N 64
int main(int argc, char *argv[]) {
  MPI_Comm world, workers;
  MPI_Group world_group, worker_group;
  int myid, nprocs;
  int server, n = -1, ranks[1];
  MPI_Init(&argc, &argv);
  world = MPI_COMM_WORLD;
  MPI_Comm_rank(world, &myid);
  MPI_Comm_size(world, &nprocs);
  server = nprocs-1;
  MPI_Comm_group(world, &world_group);
  ranks[0] = server;
  MPI_Group_excl(world_group, 1, ranks, &worker_group);
  MPI_Comm_create(world, worker_group, &workers);
  MPI_Group_free(&worker_group);
  if (myid != server) {
    MPI_Allreduce(&myid, &n, 1, MPI_INT, MPI_SUM, workers);
    MPI_Comm_free(&workers);
  }
  printf("process %2d: n = %6d\n", myid, n);
  MPI_Finalize();
  return 0;
}
