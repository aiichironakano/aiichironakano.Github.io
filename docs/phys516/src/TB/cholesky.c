#include <stdio.h>

double **dmatrix(int, int, int, int);
double *dvector(int, int);
void choldc(double **, int, double *);

int main() {
	double **a;
	double **l;
	double *p;
	int i,j,k;
	double dummy = 0.0;
	double sum;

	a = dmatrix(1,3,1,3);
	l = dmatrix(1,3,1,3);
	p = dvector(1,3);

	a[1][1] = 1.0;
	a[1][2] = 0.2;
	a[1][3] = 0.1;
	a[2][2] = 1.0;
	a[2][3] = 0.3;
	a[3][3] = 1.0;

	choldc(a,3,p);

	for (i=1; i<=3; i++) {
		for (j=1; j<i; j++)
			l[i][j] = a[i][j];
		l[i][i] = p[i];
	}

	printf("A\n");
	for (i=1; i<=3; i++) {
		for (j=1; j<i; j++)
			printf("%le\t",a[j][i]);
		printf("%le\t",a[i][i]);
		for (j=i+1; j<=3; j++)
			printf("%le\t",a[i][j]);
		printf("\n");
	}

	printf("\nL\n");
	for (i=1; i<=3; i++) {
		for (j=1; j<i; j++)
			printf("%le\t",a[i][j]);
		printf("%le\n",p[i]);
	}

	printf("\nL•Lt\n");
	for (i=1; i<=3; i++) {
		for (j=1; j<=3; j++) {
			for (k=1,sum=0.0; k<=3; k++)
				sum += l[i][k]*l[j][k];
			printf("%le\t",sum);
		}
		printf("\n");
	}

	return 0;
}
