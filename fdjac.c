#include <math.h>
#include "nrutil.h"
#define EPS 1.0e-7

void fdjac(n,x,fvec,df,vecfunc)
double **df,fvec[],x[];
int n;
void (*vecfunc)();
{
	int i,j;
	double h,temp,*f;

	f=dvector(1,n);
	for (j=1;j<=n;j++) {
		temp=x[j];
		h=EPS*fabs(temp);
		if (h == 0.0) h=EPS;
		x[j]=temp+h;
		h=x[j]-temp;
		(*vecfunc)(n,x,f);
		x[j]=temp;
		for (i=1;i<=n;i++) df[i][j]=(f[i]-fvec[i])/h;
	}
	free_dvector(f,1,n);
}
#undef EPS
