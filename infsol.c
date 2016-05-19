
#include <math.h>
#include <stdio.h>

/* find inflow solution for grmhd code,
 * generate initial conditions 
 * assumes Kerr-Schild coordinates */

double a ;
double GM = 1. ;
double gam = 5./3. ;

double DEL = 2.e-5 ;    /* optimal choice for centered difference 
                           numerical derivatives */

double FM = -1. ;	/* "rest mass" flux */
double FW ;		/* sqrt(g) (u^r b^\phi - u^\phi b^r) */
double FB ;		/* sqrt(g) b^r */
double K ;		/* entropy constant */

double r_isco,omega_isco,uphi_isco ;

int main(int argc, char *argv[])
{
	double efl(double r, double ur, double FL) ;
	void isco_calc() ;
	void find_eigenvalues() ;

	fprintf(stderr,"hello, world\n") ;

	/* set global parameters of flow */
	a = 0. ;
	FB = 1.e-7 ;
	K = 0.00001 ;

	find_eigenvalues() ;

	return(1) ;
}

void find_eigenvalues() 
{
	void isco_calc() ;
	double x[5] ;
	int nstep,i,check ;
	double FBi,FBf,af,ai,dFB,da ;
	void vecfunc(int n, double *x, double *fvec) ;
	void newt(double *x,int n, int *check, void (*vecfunc)() ) ;

	/* initial guesses on eigenvalues */
	x[0] = 5.8 ;	/* slow point radius */
	x[1] = 5.0 ;	/* fast point radius */

	x[2] = -0.001 ;	/* slow point velocity */
	x[3] = -0.004 ;	/* fast point velocity */

	x[4] = -1.9 ;	/* angular momentum flux */

	isco_calc() ;

	fprintf(stderr,"**newt called in eigen\n") ;
	newt(x-1, 5, &check, vecfunc) ;
	fprintf(stderr,"%d %g %g %g %g %g\n",i,x[0],x[1],x[2],x[3],x[4]) ;
	fprintf(stderr,"**newt done in eigen\n") ;


#if 0
	nstep = 100 ;
	FBi = 4. ;
	FBf = FB ;
	ai = 0. ;
	af = a ;
	dFB = (FBf - FBi)/nstep ;
	da = (af - ai)/nstep ;

        for(i = 0 ; i <= 3*nstep ; i++) {
		/* this is a path from _i to _f that works */
		if(i < nstep)
			a = 0.5*i/nstep ;
		else if(i < 2*nstep)
			a = 0.5 ;
		else
			a = 0.5 + (af - 0.5)*(i - 2*nstep)/nstep ;

		if(i < nstep)
			FB = 4. ;
		else if(i < 2*nstep)
			FB = 4. + (FBf - 4.)*(i - nstep)/nstep ;
		else
			FB = FBf ;

		isco_calc() ;
		fprintf(stderr,"**newt called in eigen\n") ;
		newt(x-1, 5, &check, vecfunc) ;
		fprintf(stderr,"%d %g %g %g %g %g\n",i,x[0],x[1],x[2],x[3],x[4]) ;
		fprintf(stderr,"**newt done in eigen\n") ;
	}
#endif

}

#if 0
/** once eigenvalues are found, trace out solution **/
taka_flow(double r, double *prim) 
{
	static int firstc = 1 ;

	if(firstc) {
		firstc = 0 ;
	}

}
#endif

/* set r_isco, omega_isco, uphi_isco */
void isco_calc()
{
        double Z1,Z2 ;
	double gdet_isco, br_isco ;
	extern double a ;
	extern double r_isco,omega_isco,uphi_isco ;

        /* first find rmso */
        Z1 = 1. + pow(1. - a*a,1./3.)*(pow(1. + a,1./3.) +
                pow(1. - a, 1./3.)) ;
        Z2 = sqrt(3.*a*a + Z1*Z1) ;

        r_isco = 3. + Z2 - sqrt((3. - Z1)*(3. + Z1 + 2.*Z2)) ;

        /* now calculate W_isco for prograde orbit */
        omega_isco = 1./(pow(r_isco, 1.5) + a) ;

	/* now calculate uphi_isco */
	uphi_isco = 1./sqrt( 
			1.  - 2./r_isco + 4.*a*omega_isco/r_isco
			- (a*a + 2.*a*a/r_isco + r_isco*r_isco)
			*omega_isco*omega_isco 
			) ;

	gdet_isco = r_isco*r_isco ;
	br_isco = FB/gdet_isco ;

	/* tell us about it */
	fprintf(stderr,"ISCO r, omega, uphi, gdet, br: %g %g %g %g %g\n",
			r_isco,omega_isco, uphi_isco, gdet_isco,
			br_isco) ;

	/* plausible setting for FW */
	FW = -1.0*gdet_isco*br_isco*uphi_isco ;
}

void vecfunc(int n, double *x, double *fvec)
{
        double rs,rf,us,uf,FL,Ef,Es ;
        double dedr(double r, double ur, double FL) ;
        double dedur(double r, double ur, double FL) ;
        double efl(double r, double ur, double FL) ;
        
        rs = x[1] ;
        rf = x[2] ;
        us = x[3] ;
        uf = x[4] ;
        FL = x[5] ;

        /* eqtn. 1: slow velocity condition */
        fvec[1] = dedur(rs,us,FL) ;

        /* eqtn. 2: fast velocity condition */
        fvec[2] = dedur(rf,uf,FL) ;

        /* eqtn. 3: dE/dr = 0 at slow point */
        fvec[3] = dedr(rs,us,FL) ;

        /* eqtn. 4: dE/dr = 0 at fast point */
        fvec[4] = dedr(rf,uf,FL) ;

        /* eqtn. 5: E(fast) = E(slow) */
        Es = efl(rs,us,FL) ;
        Ef = efl(rf,uf,FL) ;
        fvec[5] = Ef - Es ;


	fprintf(stderr,"rs,rf,us,uf,FL: %g %g %g %g %g\n",
			x[1],
			x[2],
			x[3],
			x[4],
			x[5]) ;
	fprintf(stderr,"vecfuncs: %g %g %g %g %g\n",
			fvec[1],
			fvec[2],
			fvec[3],
			fvec[4],
			fvec[5]) ;
}

double dedr(double r, double ur, double FL) 
{
        double del ;
        double efl(double r, double ur, double FL) ;

        del = (efl(r*(1. + DEL), ur, FL) - efl(r*(1. - DEL), ur, FL))/
                (efl(r, ur, FL)*2.*DEL) ;

        return(del) ;
}

double dedur(double r, double ur, double FL) 
{
        double del ;
        double efl(double r, double ur, double FL) ;

        del = (efl(r, ur*(1. + 0.5*DEL), FL) 
                - efl(r, ur*(1. - 0.5*DEL), FL))/
                (efl(r, ur, FL)*DEL) ;

        return(del) ;
}

double GDET,R,UT,UR,BR,RHO,U ;
double TARGET_FL ;

double efl(double r, double ur, double FL) 
{
	double x[2] ;
	int check ;
	void vecfunc2() ;
	double energy_flux_calc(double uphi, double bphi) ;
	double energy_flux ;
	void newt2(double *x,int n, int *check, void (*vecfunc)() ) ;

	/* set the easy stuff */
	R = r ;
	GDET = R*R ;
	UR = ur ;
	RHO = FM/(GDET*UR) ;
	U = K*pow(RHO,gam)/(gam - 1.) ;
	BR = FB/GDET ;

	TARGET_FL = FL ;

	/* first: find uphi, bphi given FL, FW */
	x[0] = 0.4 ; 		/* guess for uphi */
	x[1] = -BR/R ;      	/* guess for bphi */

	newt2(x-1,2,&check,vecfunc2) ;

	if(check) {
		fprintf(stderr,"newt2 failed.\n") ;
		exit(1) ;
	}

	/* now return energy flux */
	energy_flux = energy_flux_calc(x[0],x[1]) ;
	/*
	*/
	return( energy_flux ) ;

}

void vecfunc2(int n, double *x, double *vec)
{
	double up,bp ;
	double ang_mom_flux_calc(double up, double bp), utcalc(double up) ;
	extern double GDET,UT,UR,BR,FW,TARGET_FL ;

	up = x[1] ;
	bp = x[2] ;

	UT = utcalc(up) ;

	/* FL */
	vec[1] = ang_mom_flux_calc(up,bp) - TARGET_FL ;
	/* FW */
	vec[2] = (-BR*up/UT + bp*UR/UT)*GDET - FW ;

	/*
	fprintf(stderr,"vecfunc2: uphi,bphi,vec1,vec2 | %g %g %g %g\n",up,bp,vec[1],vec[2]) ;
	*/
}

double ang_mom_flux_calc(double UP, double BP)
{
	double p,upd,bt,br,bp,b2,bpd ;
	extern double a,R,U,RHO,UT,UR,BR ;

	p = (gam - 1.)*U ;
	upd = UT*(-2.*a/R) + UR*(-a*(1.+2./R)) + UP*(a*a*(1. + 2./R) + R*R) ;

	bt = BR*(-a*UP + UR + (2./R)*(-a*UP + UR + UT))
		+ BP*(a*a*UP + R*R*UP - a*UR + (2./R)*(a*a*UP - a*UR - a*UT)) ;
	br = (BR + bt*UR)/UT ;
	bp = (BP + bt*UP)/UT ;
	b2 = bt*(bt*(-1. + 2./R) - 2.*a*bp/R + 2.*br/R)
	   + br*(-(a*bp*(1. + 2./R)) + br*(1. + 2./R) + 2.*bt/R)
	   + bp*(-(a*br*(1. + 2./R)) - 2.*a*bt/R + bp*(a*a*(1. + 2./R) + R*R)) ;
	bpd = bt*(-2.*a/R) + br*(-a*(1.+2./R)) + bp*(a*a*(1.+2./R) + R*R) ;

	return( ((RHO + p + U + b2)*UR*upd - br*bpd)*GDET ) ;
}

double energy_flux_calc(double UP, double BP)
{
	double p,utd,bt,br,bp,b2,btd ;
	extern double a,R,GDET,U,RHO,UT,UR,BR ;

	p = (gam - 1.)*U ;
	utd = -2.*a*UP/R + 2.*UR/R + (-1. + 2./R)*UT ;

	bt = BR*(-a*UP + UR + (2./R)*(-a*UP + UR + UT))
		+ BP*(a*a*UP + R*R*UP - a*UR + (2./R)*(a*a*UP - a*UR - a*UT)) ;
	br = (BR + bt*UR)/UT ;
	bp = (BP + bt*UP)/UT ;
	b2 = bt*(bt*(-1. + 2./R) - 2.*a*bp/R + 2.*br/R)
	   + br*(-(a*bp*(1. + 2./R)) + br*(1. + 2./R) + 2.*bt/R)
	   + bp*(-(a*br*(1. + 2./R)) - 2.*a*bt/R + bp*(a*a*(1. + 2./R) + R*R)) ;
	btd = bt*(-1. + 2./R) - 2.*a*bp/R + 2.*br/R ;

	return( ((RHO + p + U + b2)*UR*utd - br*btd)*GDET ) ;
}

double utcalc(double up) 
{
	extern double R,UR,a ;
	double ut ;

	ut = (-2.*a*up + 2.*UR + R*sqrt(
			1. - 2./R + a*a*up*up
			- 2.*R*up*up
			+ R*R*up*up
			- 2.*a*up*UR
			+ UR*UR))/(R - 2.)  ;
	fprintf(stderr,"
	return( ut ) ;
}

