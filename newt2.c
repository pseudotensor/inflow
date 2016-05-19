

/*
 *
 * newt22: efficiently solve set of 2 nonlinear equations
 * using newt2on-raphson.
 *
 * cfg 12-14-01
 *
 */

#include <math.h>
#include <stdio.h>
#include "nrutil.h"

#define MAXITS 	200
#define TOLF 	1.0e-14
#define TOLMIN 	1.0e-14
#define TOLX 	1.0e-14
#define STPMX 	1.0

int nn2;
double *fvec2 ;
void (*nrfuncv2)();

#define FREERETURN {free_dvector(fvec2,1,n);free_dvector(xold,1,n);\
	free_dvector(p,1,n);free_dvector(g,1,n);free_dmatrix(fjac,1,n,1,n);\
	free_ivector(indx,1,n);return;}

void newt2(x,n,check,vecfunc2)
double x[];
int *check,n;
void (*vecfunc2)();
{
	double fmin2();
	void fdjac(),lnsrch(),lubksb(),ludcmp();
	int i,its,j,*indx;
	double d,den,f,fold,stpmax,sum,temp,test,**fjac,*g,*p,*xold;

	indx=ivector(1,n);
	fjac=dmatrix(1,n,1,n);

	g    = dvector(1,n);
	p    = dvector(1,n);
	xold = dvector(1,n);
	fvec2 = dvector(1,n);

	nn2=n;
	nrfuncv2=vecfunc2;
	f=fmin2(x);

	for(x[1]=-20.;x[1]<20.;x[1] += 0.2)
	for(x[2]=-20.;x[2]<20.;x[2] += 0.2) {
		fprintf(stderr,"%g %g %g\n",x[1],x[2],fmin2(x)) ;
	}
	exit(0) ;

	test=0.0;
	for (i=1;i<=n;i++)
		if (fabs(fvec2[i]) > test) test=fabs(fvec2[i]);
	if (test < 0.01*TOLF) {
		*check=0;
		FREERETURN
	}
	for (sum=0.0,i=1;i<=n;i++) sum += SQR(x[i]);
	stpmax=STPMX*DMAX(sqrt(sum),(double)n);
	for (its=1;its<=MAXITS;its++) {
		fprintf(stderr,"its2: %d\n",its) ;

		/*
		for(i=1;i<=n;i++) fprintf(stderr,"%g ",x[i]) ;
		fprintf(stderr,"\n") ;
		*/
		
		fdjac(n,x,fvec2,fjac,vecfunc2);
		for (i=1;i<=n;i++) {
			for (sum=0.0,j=1;j<=n;j++) sum += fjac[j][i]*fvec2[j];
			g[i]=sum;
		}
		for (i=1;i<=n;i++) xold[i]=x[i];
		fold=f;
		for (i=1;i<=n;i++) p[i] = -fvec2[i];
		ludcmp(fjac,n,indx,&d);
		lubksb(fjac,n,indx,p);
		lnsrch(n,xold,fold,g,p,x,&f,stpmax,check,fmin2);
		test=0.0;
		for (i=1;i<=n;i++) {
			if (fabs(fvec2[i]) > test) test=fabs(fvec2[i]);
			fprintf(stderr,"i,fv: %d %g\n",i,fvec2[i]) ;
			/*
			*/
		}
		if (test < TOLF) {
			/*
			for (i=1;i<=n;i++) {
				fprintf(stderr,">i,fv: %d %g\n",i,fvec2[i]) ;
			}
			*/
			*check=0;
			fprintf(stderr,"return a\n") ;

			FREERETURN
		}
		if (*check) {
			test=0.0;
			den=DMAX(f,0.5*n);
			fprintf(stderr,"den: %g\n",den) ;
			for (i=1;i<=n;i++) {
				temp=fabs(g[i])*DMAX(fabs(x[i]),1.0)/den;
				fprintf(stderr,"%d %g %g %g\n",i,g[i],x[i],temp) ;
				if (temp > test) test=temp;
			}
			fprintf(stderr,"%g\n",test) ;
			*check=(test < TOLMIN ? 1 : 0);

			fprintf(stderr,"return b *check %d\n",*check) ;
			FREERETURN
		}
		test=0.0;
		for (i=1;i<=n;i++) {
			temp=(fabs(x[i]-xold[i]))/DMAX(fabs(x[i]),1.0);
			if (temp > test) test=temp;
		}
		if (test < TOLX) {
			fprintf(stderr,"return c\n") ;
			FREERETURN
		}
	}
	nrerror("MAXITS exceeded in newt2");

}

#undef MAXITS
#undef TOLF
#undef TOLMIN
#undef TOLX
#undef STPMX
#undef FREERETURN

double fmin2(x)
double x[];
{
	int i;
	double sum;

	(*nrfuncv2)(nn2,x,fvec2);
	for (sum=0.0,i=1;i<=nn2;i++) sum += SQR(fvec2[i]);
	return 0.5*sum;
}

