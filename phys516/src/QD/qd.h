/*******************************************************************************
File qd.h is a header file for program qd.c.
*******************************************************************************/
#define ND 2     /* Number of spatial dimensions (x & y) = 2 */
#define NX 512   /* Number of mesh points in the x direction */
#define NY 64    /*                             y            */
#define PI 3.141592653589793

/* Function prototypes ********************************************************/
void init_param();
void init_prop();
void init_wavefn();
void single_step();
void pot_prop();
void kin_prop(int, int);
void periodic_bc();
void calc_energy();

/* Input parameters ***********************************************************/
double LX,LY;    /* Simulation box lengths in the x & y directions */
double DT;       /* Time discretization unit */
int NSTEP;       /* Number of simulation steps */
int NECAL;       /* Interval to calculate energies */
double X0,S0,E0; /* Center-of-mass, spread & energy of initial wave packet */
double BH,BW;    /* Barrier height & width */
double EH;       /* Edge potential height */

/* Arrays **********************************************************************
psi[NX+2][NY+2][2]:  psi[i][j][0|1] is the real|imaginary part of the wave
                     function on mesh point (i,j) in the xy plane
wrk[NX+2][NY+2][2]:  Work array for a wave function
al[ND][2][2]:        al[0|1][0|1][] is the x|y-direction half|full-step
                     diagonal kinetic propagator (complex)
bux|y[2][NX|Y+2][2]: bux|y[0|1][i][] is the x|y-direction half|full-step
                     upper off-diagonal kinetic propagator on mesh i (complex)
blx|y[2][NX|Y+2][2]: blx|y[0|1][i][] is the x|y-direction half|full-step
                     lower off-diagonal kinetic propagator on mesh i (complex)
v[NX+2][NY+2]:       v[i][j] is the potential energy at mesh point (i,j)
u[NX+2][NY+2][2]:    u[i][j][] is the potential propagator on (i,j) (complex)
*******************************************************************************/
double psi[NX+2][NY+2][2];
double wrk[NX+2][NY+2][2];
double al[ND][2][2];
double bux[2][NX+2][2],blx[2][NX+2][2],buy[2][NY+2][2],bly[2][NY+2][2];
double v[NX+2][NY+2];
double u[NX+2][NY+2][2];

/* Variables *******************************************************************
dx,dy = Mesh spacing in the x & y directions
ekin  = Kinetic energy
epot  = Potential energy
etot  = Total energy
*******************************************************************************/
double dx,dy;
double ekin,epot,etot;
