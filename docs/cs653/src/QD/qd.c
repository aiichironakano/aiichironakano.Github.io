/*******************************************************************************
Quantum dynamics (QD) simulation of an electron.

USAGE

%cc -o qd qd.c -lm
%qd < qd.in (see qd.h for the input-file format)
*******************************************************************************/
#include <stdio.h>
#include <math.h>
#include "qd.h"

main(int argc, char **argv) {
	int step;

	init_param();
	init_prop();
	init_wavefn();

	for (step=1; step<=NSTEP; step++) {
		single_step();
		if (step%NECAL==0) {
			calc_energy();
			printf("%le %le %le %le\n",DT*step,ekin,epot,etot);
		}
	}
}

/*----------------------------------------------------------------------------*/
void init_param() {
/*------------------------------------------------------------------------------
	Initializes parameters.
------------------------------------------------------------------------------*/
	/* Read control parameters */
	scanf("%le%le",&LX,&LY);
	scanf("%le",&DT);
	scanf("%d",&NSTEP);
	scanf("%d",&NECAL);
	scanf("%le%le%le",&X0,&S0,&E0);
	scanf("%le%le",&BH,&BW);
	scanf("%le",&EH);

	/* Calculate the mesh size */
	dx = LX/NX;
	dy = LY/NY;
}

/*----------------------------------------------------------------------------*/
void init_prop() {
/*------------------------------------------------------------------------------
	Initializes the kinetic & potential propagators.
------------------------------------------------------------------------------*/
	int dir,lim,stp,s,i,up,lw;
	int sx,sy;
	double a,exp_p[2],ep[2],em[2];
	double x,y;

	/* Set up kinetic propagators */
	for (dir=0; dir<ND; dir++) { /* Loop over x & y directions */
		if (dir==0) {
			lim=NX;
			a=0.5/(dx*dx);
		}
		else if (dir==1) {
			lim=NY;
			a=0.5/(dy*dy);
		}

		for (stp=0; stp<2; stp++) { /* Loop over half & full steps */
			exp_p[0]=cos(-(stp+1)*DT*a);
			exp_p[1]=sin(-(stp+1)*DT*a);
			ep[0]=0.5*(1.0+exp_p[0]);
			ep[1]=0.5*exp_p[1];
			em[0]=0.5*(1.0-exp_p[0]);
			em[1]= -0.5*exp_p[1];

			/* Diagonal propagator */
			for (s=0; s<2; s++) al[dir][stp][s]=ep[s];

			/* upper & lower subdiagonal propagators */
			for (i=1; i<=lim; i++) { /* Loop over mesh points */
				if (stp==0) { /* Half-step */
					up=i%2;     /* Odd mesh point has upper off-diagonal */
					lw=(i+1)%2; /* Even               lower              */
				}
				else { /* Full step */
					up=(i+1)%2; /* Even mesh point has upper off-diagonal */
					lw=i%2;     /* Odd                 lower              */
				}
				if (dir==0)
					for (s=0; s<2; s++) {
						bux[stp][i][s]=up*em[s];
						blx[stp][i][s]=lw*em[s];
					}
				else if (dir==1)
					for (s=0; s<2; s++) {
						buy[stp][i][s]=up*em[s];
						bly[stp][i][s]=lw*em[s];
					}
			} /* Endfor mesh points */
		} /* Endfor half & full steps */
	} /* Endfor x & y directions */

	/* Set up potential propagator */
	for (sx=1; sx<=NX; sx++) {
		x=dx*sx;
		for (sy=1; sy<=NY; sy++) {
			y=dy*sy;
			/* Construct the edge potential */
			if (sx==1 || sx==NX || sy==1 || sy==NY)
				v[sx][sy]=EH;
			/* Construct the barrier potential */
			else if (0.5*(LX-BW)<x && x<0.5*(LX+BW))
				v[sx][sy]=BH;
			else
				v[sx][sy]=0.0;
			/* Half-step potential propagator */
			u[sx][sy][0]=cos(-0.5*DT*v[sx][sy]);
			u[sx][sy][1]=sin(-0.5*DT*v[sx][sy]);
		}
	}
}

/*----------------------------------------------------------------------------*/
void init_wavefn() {
/*------------------------------------------------------------------------------
	Initializes the wave function as a Gaussian wave packet.
------------------------------------------------------------------------------*/
	int sx,sy,s;
	double x,y,gauss,psisq,norm_fac;

	/* Calculate the value of the wave function mesh point-by-point */
	for (sx=1; sx<=NX; sx++) {
		x = dx*sx-X0;
		for (sy=1; sy<=NY; sy++) {
			if (sy==1 || sy==NY)
				gauss = 0.0;
			else {
				y = dy*sy;
				gauss = exp(-0.25*x*x/(S0*S0))*sin(PI*(y-dy)/(LY-2.0*dy));
			}
			psi[sx][sy][0] = gauss*cos(sqrt(2.0*E0)*x);
			psi[sx][sy][1] = gauss*sin(sqrt(2.0*E0)*x);
		}
	}

	/* Normalize the wave function */
	psisq=0.0;
	for (sx=1; sx<=NX; sx++)
		for (sy=1; sy<=NY; sy++)
			for (s=0; s<=1; s++)
				psisq += psi[sx][sy][s]*psi[sx][sy][s];
	psisq *= dx*dy;
	norm_fac = 1.0/sqrt(psisq);
	for (sx=1; sx<=NX; sx++)
		for (sy=1; sy<=NY; sy++)
			for (s=0; s<=1; s++)
				psi[sx][sy][s] *= norm_fac;
}

/*----------------------------------------------------------------------------*/
void single_step() {
/*------------------------------------------------------------------------------
	Propagates the electron wave function for a unit time step, DT.
------------------------------------------------------------------------------*/
	pot_prop();    /* half step potential propagation                  */

	kin_prop(0,0); /* half step kinetic propagation in the x direction */
	kin_prop(0,1); /* full                                             */
	kin_prop(0,0); /* half                                             */

	kin_prop(1,0); /* half step kinetic propagation in the y direction */
	kin_prop(1,1); /* full                                             */
	kin_prop(1,0); /* half                                             */

	pot_prop();    /* half step potential propagation                  */
}

/*----------------------------------------------------------------------------*/
void pot_prop() {
/*------------------------------------------------------------------------------
	Potential propagator for a half time step, DT/2.
------------------------------------------------------------------------------*/
	int sx,sy;
	double wr,wi;

	for (sx=1; sx<=NX; sx++)
		for (sy=1; sy<=NY; sy++) {
			wr=u[sx][sy][0]*psi[sx][sy][0]-u[sx][sy][1]*psi[sx][sy][1];
			wi=u[sx][sy][0]*psi[sx][sy][1]+u[sx][sy][1]*psi[sx][sy][0];
			psi[sx][sy][0]=wr;
			psi[sx][sy][1]=wi;
		}
}

/*----------------------------------------------------------------------------*/
void kin_prop(int d, int t) {
/*------------------------------------------------------------------------------
	Kinetic propagation in the d (=0 for x; 1 for y) direction
	for t (=0 for DT/2--half; 1 for DT--full) step.
-------------------------------------------------------------------------------*/
	int sx,sy,s;
	double wr,wi;

	/* Apply the periodic boundary condition */
	periodic_bc();

	/* WRK|PSI holds the new|old wave function */
	for (sx=1; sx<=NX; sx++)
		for (sy=1; sy<=NY; sy++) {
			wr=al[d][t][0]*psi[sx][sy][0]-al[d][t][1]*psi[sx][sy][1];
			wi=al[d][t][0]*psi[sx][sy][1]+al[d][t][1]*psi[sx][sy][0];
			if (d==0) {
				wr+=(blx[t][sx][0]*psi[sx-1][sy][0]-blx[t][sx][1]*psi[sx-1][sy][1]);
				wi+=(blx[t][sx][0]*psi[sx-1][sy][1]+blx[t][sx][1]*psi[sx-1][sy][0]);
				wr+=(bux[t][sx][0]*psi[sx+1][sy][0]-bux[t][sx][1]*psi[sx+1][sy][1]);
				wi+=(bux[t][sx][0]*psi[sx+1][sy][1]+bux[t][sx][1]*psi[sx+1][sy][0]);
			}
			else if (d==1) {
				wr+=(bly[t][sy][0]*psi[sx][sy-1][0]-bly[t][sy][1]*psi[sx][sy-1][1]);
				wi+=(bly[t][sy][0]*psi[sx][sy-1][1]+bly[t][sy][1]*psi[sx][sy-1][0]);
				wr+=(buy[t][sy][0]*psi[sx][sy+1][0]-buy[t][sy][1]*psi[sx][sy+1][1]);
				wi+=(buy[t][sy][0]*psi[sx][sy+1][1]+buy[t][sy][1]*psi[sx][sy+1][0]);
			}
			wrk[sx][sy][0]=wr;
			wrk[sx][sy][1]=wi;
		}

	/* Copy the new wave function back to PSI */
	for (sx=1; sx<=NX; sx++)
		for (sy=1; sy<=NY; sy++)
			for (s=0; s<=1; s++)
				psi[sx][sy][s]=wrk[sx][sy][s];
}

/*----------------------------------------------------------------------------*/
void periodic_bc() {
/*------------------------------------------------------------------------------
	Applies the periodic boundary condition to wave function PSI.
------------------------------------------------------------------------------*/
	int sx,sy,s;

	/* Copy boundary wave function values in the x direction */
	for (sy=1; sy<=NY; sy++)
		for (s=0; s<=1; s++) {
			psi[0][sy][s] = psi[NX][sy][s];
			psi[NX+1][sy][s] = psi[1][sy][s];
		}

	/* Copy boundary wave function values in the y direction */
	for (sx=1; sx<=NX; sx++)
		for (s=0; s<=1; s++) {
			psi[sx][0][s] = psi[sx][NY][s];
			psi[sx][NY+1][s] = psi[sx][1][s];
		}
}

/*----------------------------------------------------------------------------*/
void calc_energy() {
/*------------------------------------------------------------------------------
	Calculates the kinetic, potential & total energies.
------------------------------------------------------------------------------*/
	int sx,sy,s;
	double a,bx,by;

	/* Apply the periodic boundary condition */
	periodic_bc();

	/* Tridiagonal kinetic-energy operators */
	a=1.0/(dx*dx)+1.0/(dy*dy);
	bx= -0.5/(dx*dx);
	by= -0.5/(dy*dy);

	/* |WRK> = (-1/2)Laplacian|PSI> */
	for (sx=1; sx<=NX; sx++)
		for (sy=1; sy<=NY; sy++)
			for (s=0; s<=1; s++)
				wrk[sx][sy][s] = a*psi[sx][sy][s]
				               + bx*(psi[sx-1][sy][s]+psi[sx+1][sy][s])
				               + by*(psi[sx][sy-1][s]+psi[sx][sy+1][s]);

	/* Kinetic energy = <PSI|(-1/2)Laplacian|PSI> = <PSI|WRK> */
	ekin = 0.0;
	for (sx=1; sx<=NX; sx++)
		for (sy=1; sy<=NY; sy++)
			ekin += (psi[sx][sy][0]*wrk[sx][sy][0]+psi[sx][sy][1]*wrk[sx][sy][1]);
	ekin *= dx*dy;

	/* Potential energy */
	epot = 0.0;
	for (sx=1; sx<=NX; sx++)
		for (sy=1; sy<=NY; sy++)
			epot += v[sx][sy]*(psi[sx][sy][0]*psi[sx][sy][0]+psi[sx][sy][1]*psi[sx][sy][1]);
	epot *= dx*dy;

	/* Total energy */
	etot = ekin+epot;
}
