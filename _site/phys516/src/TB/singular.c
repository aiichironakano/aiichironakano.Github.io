#include <stdio.h>

double **dmatrix(int, int, int, int);
double *dvector(int, int);
void svdcmp(double **, int, int, double *, double **);

int main() {
	double **a;
	double *w;
	double **u,**v;
	int i,j,k;
	double t;
	double t1[3],t2[3];

	a = dmatrix(1,3,1,3);
	u = dmatrix(1,3,1,3);
	w = dvector(1,3);
	v = dmatrix(1,3,1,3);

	a[1][1] = 3.0;
	a[1][2] = 2.0;
	a[1][3] = 1.0;
	a[2][2] = 2.0;
	a[2][3] = 0.0;
	a[3][3] = 1.0;
	a[2][1] = a[1][2];
	a[3][1] = a[1][3];
	a[3][2] = a[2][3];

	for (i=1; i<=3; i++) {
		for (j=1; j<=3; j++)
			u[i][j] = a[i][j];
		printf("%le %le %le\n",u[i][1],u[i][2],u[i][3]);
	}

	svdcmp(u,3,3,w,v);

	/* Sort the singular values in descending order */
	for (i=1; i<3; i++) {
		for (j=i+1; j<=3; j++) {
			if (w[i]<w[j]) {
				t = w[i];
				w[i] = w[j];
				w[j] = t;
				for (k=1; k<=3; k++) {
					t1[k] = u[k][i];
					t2[k] = v[k][i];
				}
				for (k=1; k<=3; k++) {
					u[k][i] = u[k][j];
					v[k][i] = v[k][j];
				}
				for (k=1; k<=3; k++) {
					u[k][j] = t1[k];
					v[k][j] = t2[k];
				}
			}
		}
	}

	for (i=1; i<=3; i++) {
		printf("Sigular value %d = %le\n",i,w[i]);
		printf("        vector   = %le %le %le\n",u[1][i],u[2][i],u[3][i]);
	}

	return 0;
}
