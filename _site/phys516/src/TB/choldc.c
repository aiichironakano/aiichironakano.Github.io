/*******************************************************************************
Cholesky decomposition program, svdcmp, from "Numerical Recipes in C"
(Cambridge Univ. Press) by W.H. Press, S.A. Teukolsky, W.T. Vetterling,
and B.P. Flannery
*******************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define NR_END 1
#define FREE_ARG char*
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
static double dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
(dmaxarg1) : (dmaxarg2))
static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
(iminarg1) : (iminarg2))

double **dmatrix(int nrl, int nrh, int ncl, int nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	int i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;
	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	m += NR_END;
	m -= nrl;
	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	/* return pointer to array of pointers to rows */
	return m;
}

double *dvector(int nl, int nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;
	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	return v-nl+NR_END;
}

void free_dvector(double *v, int nl, int nh)
/* free a double vector allocated with dvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

double pythag(double a, double b)
/* compute (a2 + b2)^1/2 without destructive underflow or overflow */
{
	double absa,absb;
	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb) return absa*sqrt(1.0+(absb/absa)*(absb/absa));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+(absa/absb)*(absa/absb)));
}

/******************************************************************************/
void choldc(double **a, int n, double p[])
/*******************************************************************************
Given a positive-definite symmetric matrix a[1..n][1..n], this routine
constructs its Cholesky decomposition, A = L · LT. On input, only the upper
triangle of a need be given; it is not modified. The Cholesky factor L is
returned in the lower triangle of a, except for its diagonal elements which are
returned in p[1..n].
*******************************************************************************/
{
	int i,j,k;
	double sum;

	for (i=1; i<=n; i++) {
		for (j=i; j<=n; j++) {
			for (sum=a[i][j],k=i-1; k>=1; k--)
				sum -= a[i][k]*a[j][k];
			if (i == j) {
				if (sum <= 0.0) {
					printf("choldc failed\n"); // a is not positive definite.
					exit(1);
				}
				p[i] = sqrt(sum);
			}
			else
				a[j][i] = sum/p[i];
		}
	}
}