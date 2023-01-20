/*----------------------------------------------------------------------
Program pdf0.c computes a pair distribution function for n atoms
given the 3D coordinates of the atoms.
----------------------------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

#define NHBIN 2000  // Histogram size

float al[3];        // Simulation box lengths
int n;              // Number of atoms
float *r;           // Atomic position array
FILE *fp;

float SignR(float v,float x) {if (x > 0) return v; else return -v;}

/*--------------------------------------------------------------------*/
void histogram() {
/*----------------------------------------------------------------------
Constructs a histogram NHIS for atomic-pair distribution.
----------------------------------------------------------------------*/
  float alth[3];
  float* nhis;  // Histogram array
  float rhmax,drh,rij,dr,density,gr;
  int a,ih,i,j;  

  /* Half the simulation box size */
  for (a=0; a<3; a++) alth[a] = 0.5*al[a];
  /* Max. pair distance RHMAX & histogram bin size DRH */
  rhmax = sqrt(alth[0]*alth[0]+alth[1]*alth[1]+alth[2]*alth[2]);
  drh = rhmax/NHBIN;  // Histogram bin size

  nhis = (float*)malloc(sizeof(float)*NHBIN);
  for (ih=0; ih<NHBIN; ih++) nhis[ih] = 0.0; // Reset the histogram

  for (i=0; i<n-1; i++) {
    for (j=i+1; j<n; j++) {
      rij = 0.0;
      for (a=0; a<3; a++) {
        dr = r[3*i+a]-r[3*j+a];
        /* Periodic boundary condition */
        dr = dr-SignR(alth[a],dr-alth[a])-SignR(alth[a],dr+alth[a]);
        rij += dr*dr;
      }
      rij = sqrt(rij); /* Pair distance */
      ih = rij/drh;
      nhis[ih] += 1.0; /* Entry to the histogram */
    }  // End for j
  }  // Endo for i

  density = n/(al[0]*al[1]*al[2]);
  /* Print out the histogram */
  fp = fopen("pdf.d","w");
  for (ih=0; ih<NHBIN; ih++) {
    gr = nhis[ih]/(2*M_PI*pow((ih+0.5)*drh,2)*drh*density*n);
    fprintf(fp,"%e %e\n",(ih+0.5)*drh,gr);
  }
  fclose(fp);
  free(nhis);
}

/*--------------------------------------------------------------------*/
int main() {
/*--------------------------------------------------------------------*/
  int i;
  float cpu1,cpu2;

  /* Read the atomic position data */
  fp = fopen("pos.d","r");
  fscanf(fp,"%f %f %f",&(al[0]),&(al[1]),&(al[2]));
  fscanf(fp,"%d",&n);
  r = (float*)malloc(sizeof(float)*3*n);
  for (i=0; i<n; i++)
    fscanf(fp,"%f %f %f",&(r[3*i]),&(r[3*i+1]),&(r[3*i+2]));
  fclose(fp);

  /* Compute the histogram */
  cpu1 = ((float) clock())/CLOCKS_PER_SEC;
  histogram();  
  cpu2 = ((float) clock())/CLOCKS_PER_SEC;
  printf("Execution time (s) = %le\n",cpu2-cpu1);

  free(r);
  return 0;
}
