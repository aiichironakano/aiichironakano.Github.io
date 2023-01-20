/*******************************************************************************
File md.h is an include file for program md.c.
*******************************************************************************/
#define NMAX 100000  /* Maximum number of atoms which can be simulated */
#define RCUT 2.5     /* Potential cut-off length */
#define PI 3.141592653589793
/* Constants for the random number generator */
#define D2P31M 2147483647.0
#define DMUL 16807.0

/* Functions & function prototypes ********************************************/

double SignR(double v,double x) {if (x > 0) return v; else return -v;}
double Dmod(double a, double b) {
	int n;
	n = (int) (a/b);
	return (a - b*n);
}
double RandR(double *seed) {
	*seed = Dmod(*seed*DMUL,D2P31M);
	return (*seed/D2P31M);
}
void RandVec3(double *p, double *seed) {
	double x,y,s = 2.0;
	while (s > 1.0) {
		x = 2.0*RandR(seed) - 1.0; y = 2.0*RandR(seed) - 1.0; s = x*x + y*y;
	}
	p[2] = 1.0 - 2.0*s; s = 2.0*sqrt(1.0 - s); p[0] = s*x; p[1] = s*y;
}
void InitParams();
void InitConf();
void ComputeAccel();
void SingleStep();
void HalfKick();
void ApplyBoundaryCond();
void EvalProps();

/* Input parameters (read from an input file in this order) *******************/

int InitUcell[3];   /* Number of unit cells */
double Density;     /* Number density of atoms (in reduced unit) */
double InitTemp;    /* Starting temperature (in reduced unit) */
double DeltaT;      /* Size of a time step (in reduced unit) */
int StepLimit;      /* Number of time steps to be simulated */
int StepAvg;        /* Reporting interval for statistical data */

/* Constants ******************************************************************/

double Region[3];  /* MD box lengths */
double RegionH[3]; /* Half the box lengths */
double DeltaTH;    /* Half the time step */
double Uc, Duc;    /* Potential cut-off parameters */

/* Variables ******************************************************************/

int nAtom;            /* Number of atoms */
double r[NMAX][3];    /* r[i][0|1|2] is the x|y|z coordinate of atom i */
double rv[NMAX][3];   /* Atomic velocities */
double ra[NMAX][3];   /* Acceleration on atoms */
double kinEnergy;     /* Kinetic energy */
double potEnergy;     /* Potential energy */
double totEnergy;     /* Total energy */
double temperature;   /* Current temperature */
int stepCount;        /* Current time step */
/******************************************************************************/

