/*******************************************************************************
File fmm2d.h is an include file for program fmm2d.c.
*******************************************************************************/
#define BOX 1.0          /* Cubic simulation-box size */
#define Npar  1000       /* # of charged particles */
#define L 4              /* Max level of refinement = quadtree height */
#define P 6              /* Order of multipole & local expansions */
#define Max_par 1000000  /* Array size for particles */
#define Max_cell 1000000 /* Array size for quadtree cells */
#define Max_level 20     /* Array size for quadtree levels (>= L+1) */
#define Max_term 20      /* Array size for multipole|local terms (>= P) */
#define EMPTY -1         /* NULL pointer in linked lists */

/* Variables for timing measurement *******************************************/

double tfmm;             /* Wall-clock time for FMM */
double tdirect;          /* Wall-clock time for direct calculation */

/* Functions & function prototypes ********************************************/

void initialize();
void mp_leaf();
void upward();
void downward();
void nn_direct();
void all_direct();

/* Complex-arithmetic functions */
void cadd(double s,double *a,double t,double *b,double *c) {  /* C = sA+tB */
  c[0] = s*a[0]+t*b[0]; c[1] = s*a[1]+t*b[1];}
void smul(double *a,double s,double *c) {  /* C = sA */
  c[0] = s*a[0]; c[1] = s*a[1];}
void cmul(double *a,double *b,double *c) {  /* C = AB */
  double w[2];
  w[0] = a[0]*b[0]-a[1]*b[1];
  w[1] = a[0]*b[1]+a[1]*b[0];
  c[0] = w[0]; c[1] = w[1];
}
void cinv(double *a,double *ai) {  /* AI = 1/A */
  double a2i;
  a2i = 1.0/(a[0]*a[0]+a[1]*a[1]);
  ai[0] =  a[0]*a2i;
  ai[1] = -a[1]*a2i;
}
void clgn(double *a,double *l) {  /* L = log(A) */
  l[0] = log(sqrt(a[0]*a[0]+a[1]*a[1])); l[1] = atan(a[1]/a[0]);}
void cini(double s,double t,double *a) {  /* A = s+it */
  a[0] = s; a[1] = t;}

/* Combinatorics functions */
int fact(int n) {
  int i,val;
  if (n == 0)
    val = 1;
  else
    for (val = 1,i=1; i<=n; i++) val *= i;
  return val;
}
int comb(int n,int k) {return fact(n)/fact(k)/fact(n-k);}

/* Variables ******************************************************************/

double z[Max_par][2]; /* z[j][0|1] is the x|y coordinate of particle j */
double q[Max_par];    /* q[j] is the charge of particle j */
double phi[Max_cell][Max_term][2];
                      /* phi[c][a][0|1] is the real|imaginary part of the
                         a-th order multipole of cell c */
double psi[Max_cell][Max_term][2];
                      /* psi[c][a][0|1] is the real|imaginary part of the
                         a-th order local-expansion term of cell c */
int c0[Max_level];    /* c0[l] is the cell-index offset at level l */
double pot[Max_par];  /* Potential at j-th particle position */
double pot_direct[Max_par]; 
                      /* Potential by all-pair direct calculation */
double eng;           /* Total electrostatic energy */
double eng_direct;    /* Electrostatic energy by all-pair direct calculation */
int head[Max_cell];   /* Headers for linked lists */
int lscl[Max_par];    /* Linked lists for particles in the leaf cells */
/******************************************************************************/
