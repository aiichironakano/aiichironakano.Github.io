/* Testing Gauss Legendre integration */
#include <stdio.h>
#include <math.h>

double *dvector(int, int);
void gauleg(double, double, double *, double *, int);

int main() {
	double *x,*w;
	double x1= -1.0,x2=1.0,sum;
	int N,i;
	printf("Input the number of quadrature points\n");
	scanf("%d",&N);
	x=dvector(1,N);
	w=dvector(1,N);
	gauleg(x1,x2,x,w,N);
	sum=0.0;
	for (i=1; i<=N; i++)
		sum += w[i]*2.0/(1.0 + x[i]*x[i]);
	printf("Integration = %f\n", sum);
}
