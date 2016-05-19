#include <math.h>
#include <stdio.h>
#include "nrutil.h"

#define MAXITS 	200
#define TOLF 	1.0e-14
#define TOLMIN 	1.0e-14
#define TOLX 	1.0e-14
#define STPMX 	1.0

int nn;
double *fvec ;
void (*nrfuncv)();

#define FREERETURN {free_dvector(fvec,1,n);free_dvector(xold,1,n);\
	free_dvector(p,1,n);free_dvector(g,1,n);free_dmatrix(fjac,1,n,1,n);\
	free_ivector(indx,1,n);return;}

void newt(x,n,check,vecfunc)
double x[];
int *check,n;
void (*vecfunc)();
{
	double fmin();
	void fdjac(),lnsrch(),lubksb(),ludcmp();
	int i,its,j,*indx;
	double d,den,f,fold,stpmax,sum,temp,test,**fjac,*g,*p,*xold;

	indx=ivector(1,n);
	fjac=dmatrix(1,n,1,n);

	g    = dvector(1,n);
	p    = dvector(1,n);
	xold = dvector(1,n);
	fvec = dvector(1,n);

	nn=n;
	nrfuncv=vecfunc;
	f=fmin(x);
	test=0.0;
	for (i=1;i<=n;i++)
		if (fabs(fvec[i]) > test) test=fabs(fvec[i]);
	if (test < 0.01*TOLF) {
		*check=0;
		FREERETURN
	}
	for (sum=0.0,i=1;i<=n;i++) sum += SQR(x[i]);
	stpmax=STPMX*FMAX(sqrt(sum),(double)n);
	for (its=1;its<=MAXITS;its++) {
		if(n > 3) fprintf(stderr,"its: %d\n",its) ;

		/*
		for(i=1;i<=n;i++) fprintf(stderr,"%g ",x[i]) ;
		fprintf(stderr,"\n") ;
		*/
		
		fdjac(n,x,fvec,fjac,vecfunc);
		for (i=1;i<=n;i++) {
			for (sum=0.0,j=1;j<=n;j++) sum += fjac[j][i]*fvec[j];
			g[i]=sum;
		}
		for (i=1;i<=n;i++) xold[i]=x[i];
		fold=f;
		for (i=1;i<=n;i++) p[i] = -fvec[i];
		ludcmp(fjac,n,indx,&d);
		lubksb(fjac,n,indx,p);
		lnsrch(n,xold,fold,g,p,x,&f,stpmax,check,fmin);
		test=0.0;
		for (i=1;i<=n;i++) {
			if (fabs(fvec[i]) > test) test=fabs(fvec[i]);
			/*
			fprintf(stderr,"i,fv: %d %g\n",i,fvec[i]) ;
			*/
		}
		if (test < TOLF) {
			/*
			for (i=1;i<=n;i++) {
				fprintf(stderr,">i,fv: %d %g\n",i,fvec[i]) ;
			}
			*/
			*check=0;
			FREERETURN
		}
		if (*check) {
			test=0.0;
			den=FMAX(f,0.5*n);
			for (i=1;i<=n;i++) {
				temp=fabs(g[i])*FMAX(fabs(x[i]),1.0)/den;
				if (temp > test) test=temp;
			}
			*check=(test < TOLMIN ? 1 : 0);

			FREERETURN
		}
		test=0.0;
		for (i=1;i<=n;i++) {
			temp=(fabs(x[i]-xold[i]))/FMAX(fabs(x[i]),1.0);
			if (temp > test) test=temp;
		}
		if (test < TOLX) FREERETURN
	}
	nrerror("MAXITS exceeded in newt");

}
#undef MAXITS
#undef TOLF
#undef TOLMIN
#undef TOLX
#undef STPMX
#undef FREERETURN
