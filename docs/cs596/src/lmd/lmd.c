/*******************************************************************************
Molecular dynamics (MD) simulation with the Lennard-Jones potential with the
linked-list cell method.

USAGE

%cc -o lmd lmd.c -lm
%lmd < lmd.in (see lmd.h for the input-file format)
*******************************************************************************/
#include <stdio.h>
#include <math.h>
#include "lmd.h"
#include <time.h>

int main(int argc, char **argv) {
	double cpu, cpu1, cpu2;
	InitParams();
	InitConf();
	ComputeAccel(); /* Computes initial accelerations */
	cpu1 = ((double) clock())/CLOCKS_PER_SEC;
	for (stepCount=1; stepCount<=StepLimit; stepCount++) {
		SingleStep();
		if (stepCount%StepAvg == 0) EvalProps();
	}
	cpu2 = ((double) clock())/CLOCKS_PER_SEC;
	cpu = cpu2-cpu1;
	printf("Execution time (s) = %le\n",cpu);
	return 0;
}

/*----------------------------------------------------------------------------*/
void InitParams() {
/*------------------------------------------------------------------------------
	Initializes parameters.
------------------------------------------------------------------------------*/
	int k;
	double rr,ri2,ri6,r1;

	/* Reads control parameters */
	scanf("%d%d%d",&InitUcell[0],&InitUcell[1],&InitUcell[2]);
	scanf("%le",&Density);
	scanf("%le",&InitTemp);
	scanf("%le",&DeltaT);
	scanf("%d",&StepLimit);
	scanf("%d",&StepAvg);

	/* Computes basic parameters */
	DeltaTH = 0.5*DeltaT;
	for (k=0; k<3; k++) {
		Region[k] = InitUcell[k]/pow(Density/4.0,1.0/3.0);
		RegionH[k] = 0.5*Region[k];
	}

	/* Compute the # of cells for linked cell lists */
	for (k=0; k<3; k++) {
		lc[k] = Region[k]/RCUT; 
		rc[k] = Region[k]/lc[k];
	}

	/* Constants for potential truncation */
	rr = RCUT*RCUT; ri2 = 1.0/rr; ri6 = ri2*ri2*ri2; r1=sqrt(rr);
	Uc = 4.0*ri6*(ri6 - 1.0);
	Duc = -48.0*ri6*(ri6 - 0.5)/r1;
}

/*----------------------------------------------------------------------------*/
void InitConf() {
/*------------------------------------------------------------------------------
	r are initialized to face-centered cubic (fcc) lattice positions.
	rv are initialized with a random velocity corresponding to Temperature.
------------------------------------------------------------------------------*/
	double c[3],gap[3],e[3],vSum[3],vMag;
	int j,n,k,nX,nY,nZ;
	double seed;
	/* FCC atoms in the original unit cell */
	double origAtom[4][3] = {{0.0, 0.0, 0.0}, {0.0, 0.5, 0.5},
	                         {0.5, 0.0, 0.5}, {0.5, 0.5, 0.0}}; 

	/* Sets up a face-centered cubic (fcc) lattice */
	for (k=0; k<3; k++) gap[k] = Region[k]/InitUcell[k];
	nAtom = 0;
	for (nZ=0; nZ<InitUcell[2]; nZ++) {
		c[2] = nZ*gap[2];
		for (nY=0; nY<InitUcell[1]; nY++) {
			c[1] = nY*gap[1];
			for (nX=0; nX<InitUcell[0]; nX++) {
				c[0] = nX*gap[0];
				for (j=0; j<4; j++) {
					for (k=0; k<3; k++)
						r[nAtom][k] = c[k] + gap[k]*origAtom[j][k];
					++nAtom;
				}
			}
		}
	}

	/* Generates random velocities */
	seed = 13597.0;
	vMag = sqrt(3*InitTemp);
	for(k=0; k<3; k++) vSum[k] = 0.0;
	for(n=0; n<nAtom; n++) {
		RandVec3(e,&seed);
		for (k=0; k<3; k++) {
			rv[n][k] = vMag*e[k];
			vSum[k] = vSum[k] + rv[n][k];
		}
	}
	/* Makes the total momentum zero */
	for (k=0; k<3; k++) vSum[k] = vSum[k]/nAtom;
	for (n=0; n<nAtom; n++) for(k=0; k<3; k++) rv[n][k] = rv[n][k] - vSum[k];
}

/*----------------------------------------------------------------------------*/
void ComputeAccel() {
/*------------------------------------------------------------------------------
	Acceleration, ra, are computed as a function of atomic coordinates, r,
	using the Lennard-Jones potential.  The sum of atomic potential energies,
	potEnergy, is also computed.
------------------------------------------------------------------------------*/
	int i,j,a,lcyz,lcxyz,mc[3],c,mc1[3],c1;
	double dr[3],rr,ri2,ri6,r1,rrCut,fcVal,f,rshift[3];

	/* Reset the potential & forces */
	for (i=0; i<nAtom; i++) for (a=0; a<3; a++) ra[i][a] = 0.0;
	potEnergy = 0.0;

  /* Make a linked-cell list, lscl--------------------------------------------*/

	lcyz = lc[1]*lc[2];
	lcxyz = lc[0]*lcyz;

	/* Reset the headers, head */
	for (c=0; c<lcxyz; c++) head[c] = EMPTY;

	/* Scan atoms to construct headers, head, & linked lists, lscl */

	for (i=0; i<nAtom; i++) {
		for (a=0; a<3; a++) mc[a] = r[i][a]/rc[a];

		/* Translate the vector cell index, mc, to a scalar cell index */
		c = mc[0]*lcyz+mc[1]*lc[2]+mc[2];

		/* Link to the previous occupant (or EMPTY if you're the 1st) */
		lscl[i] = head[c];

		/* The last one goes to the header */
		head[c] = i;
	} /* Endfor atom i */

  /* Calculate pair interaction-----------------------------------------------*/

	rrCut = RCUT*RCUT;

	/* Scan inner cells */
	for (mc[0]=0; mc[0]<lc[0]; (mc[0])++)
	for (mc[1]=0; mc[1]<lc[1]; (mc[1])++)
	for (mc[2]=0; mc[2]<lc[2]; (mc[2])++) {

		/* Calculate a scalar cell index */
		c = mc[0]*lcyz+mc[1]*lc[2]+mc[2];
		/* Skip this cell if empty */
		if (head[c] == EMPTY) continue;

		/* Scan the neighbor cells (including itself) of cell c */
		for (mc1[0]=mc[0]-1; mc1[0]<=mc[0]+1; (mc1[0])++)
		for (mc1[1]=mc[1]-1; mc1[1]<=mc[1]+1; (mc1[1])++)
		for (mc1[2]=mc[2]-1; mc1[2]<=mc[2]+1; (mc1[2])++) {
			/* Periodic boundary condition by shifting coordinates */
			for (a=0; a<3; a++) {
				if (mc1[a] < 0)
					rshift[a] = -Region[a];
				else if (mc1[a]>=lc[a])
					rshift[a] = Region[a];
				else
					rshift[a] = 0.0;
			}
			/* Calculate the scalar cell index of the neighbor cell */
			c1 = ((mc1[0]+lc[0])%lc[0])*lcyz
			    +((mc1[1]+lc[1])%lc[1])*lc[2]
			    +((mc1[2]+lc[2])%lc[2]);
			/* Skip this neighbor cell if empty */
			if (head[c1] == EMPTY) continue;

			/* Scan atom i in cell c */
			i = head[c];
			while (i != EMPTY) {

				/* Scan atom j in cell c1 */
				j = head[c1];
				while (j != EMPTY) {

					/* Avoid double counting of pairs */
					if (i < j) {
						/* Pair vector dr = r[i]-r[j] */
						for (rr=0.0, a=0; a<3; a++) {
							dr[a] = r[i][a]-(r[j][a]+rshift[a]);
							rr += dr[a]*dr[a];
						}

						/* Calculate potential & forces if rij<RCUT */
						if (rr < rrCut) {
							ri2 = 1.0/rr; ri6 = ri2*ri2*ri2; r1 = sqrt(rr);
							fcVal = 48.0*ri2*ri6*(ri6-0.5) + Duc/r1;
							for (a=0; a<3; a++) {
								f = fcVal*dr[a];
								ra[i][a] += f;
								ra[j][a] -= f;
							}
							potEnergy += 4.0*ri6*(ri6-1.0) - Uc - Duc*(r1-RCUT);
						}
					} /* Endif i<j */

					j = lscl[j];
				} /* Endwhile j not empty */

				i = lscl[i];
			} /* Endwhile i not empty */

		} /* Endfor neighbor cells, c1 */
	} /* Endfor central cell, c */
}

/*----------------------------------------------------------------------------*/
void SingleStep() {
/*------------------------------------------------------------------------------
	r & rv are propagated by DeltaT in time using the velocity-Verlet method.
------------------------------------------------------------------------------*/
	int n,k;

	HalfKick(); /* First half kick to obtain v(t+Dt/2) */
	for (n=0; n<nAtom; n++) /* Update atomic coordinates to r(t+Dt) */
		for (k=0; k<3; k++) r[n][k] = r[n][k] + DeltaT*rv[n][k];
	ApplyBoundaryCond();
	ComputeAccel(); /* Computes new accelerations, a(t+Dt) */
	HalfKick(); /* Second half kick to obtain v(t+Dt) */
}

/*----------------------------------------------------------------------------*/
void HalfKick() {
/*------------------------------------------------------------------------------
	Accelerates atomic velocities, rv, by half the time step.
------------------------------------------------------------------------------*/
	int n,k;
	for (n=0; n<nAtom; n++)
		for (k=0; k<3; k++) rv[n][k] = rv[n][k] + DeltaTH*ra[n][k];
}

/*----------------------------------------------------------------------------*/
void ApplyBoundaryCond() {
/*------------------------------------------------------------------------------
	Applies periodic boundary conditions to atomic coordinates.
------------------------------------------------------------------------------*/
	int n,k;
	for (n=0; n<nAtom; n++) 
		for (k=0; k<3; k++) 
			r[n][k] = r[n][k] - SignR(RegionH[k],r[n][k])
			                  - SignR(RegionH[k],r[n][k]-Region[k]);
}

/*----------------------------------------------------------------------------*/
void EvalProps() {
/*------------------------------------------------------------------------------
	Evaluates physical properties: kinetic, potential & total energies.
------------------------------------------------------------------------------*/
	double vv;
	int n,k;

	kinEnergy = 0.0;
	for (n=0; n<nAtom; n++) {
		vv = 0.0;
		for (k=0; k<3; k++)
			vv = vv + rv[n][k]*rv[n][k];
		kinEnergy = kinEnergy + vv;
	}
	kinEnergy *= (0.5/nAtom);
	potEnergy /= nAtom;
	totEnergy = kinEnergy + potEnergy;
	temperature = kinEnergy*2.0/3.0;

	/* Print the computed properties */
	printf("%9.6f %9.6f %9.6f %9.6f\n",
	stepCount*DeltaT,temperature,potEnergy,totEnergy);
}
