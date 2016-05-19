

/* 

inf_const.c: Find eigenvalues (energy and angular
momentum flux) for fully relativistic inflow problem.

Arguments:

Fhp: (F_{\theta\phi}) component of Maxwell tensor--
	this is a magnetization parameter.  A typical
	value is somewhere around 1.
a: spin of black hole--
	works well for a ~ 0.5, but can be massaged 
	into working for a up to 0.999.
[nstep]: the code finds a series of solutions in the
	Fhp-a plane, finding (3 x nstep) solutions along
	the way.  By default nstep is 10.  The code
	starts from guesses for the eigenvalues (see
	below) that are appropriate for a = 0, Fhp = 4.
	If values are chosen that are far from this
	a larger value of nstep should be chosen.

Output: 

FL (radial angular momentum flux) 
FE (radial energy flux) 

for specified values of the mass flux (FM), the rotation 
frequency (W), and the magnetic flux (Fhp) in the fully general
relativistic inflow problem.

Output is a line suitable for input to inf_solv, which
traces out the solution.

Notation:

In general t -> time component
	   r -> radial component
	   h -> theta component
	   p -> phi component
So Ftr is $F_{tr}$ and Fhp is $F_{\theta\phi}$.


NOTES: 

*performance depends strongly on initial guess.
*uses numerical derivatives for critical point conditions
*some values of Fhp and a will simply fail.  Values that
are known to work are Fhp = 0.5, a = 0.5

CFG 3-17-99

7-17-02: updated for use in testing grmhd code

*/

#include <stdio.h>
#include <math.h>

/** global variables **/

double FM = -1. ;	/* mass flux, normalized to 1 */
double GM = 1. ;	
double Fhp ;		/* theta-phi component of Maxwell; reduces to -r^2 Br
			   in nonrelativistic limit */
double W ;		/* a constant related to the "rotation frequency
			   of the magnetic field"; here set to rotation
			   frequency at the marginally stable orbit */
double a ;		/* black hole spin parameter, |a| < 1.
			   currently prograde orbits around prograde hole
			   are assumed. */
double rmso ;		/* radius of marginally stable orbit */
double DEL = 2.e-5 ;	/* optimal choice for centered difference 
			   numerical derivatives */

int main(int argc, char *argv[])
{
	double x[3] ;
	double Fhpi,Fhpf,dFhp,ai,af,da ;
	int check ;
	void wmsocalc() ;
	void newt(double x[],int n,int *check,void (*vecfunc)()) ;
	void vecfunc(int n, double *x, double *fvec) ;
	double efl(double r, double ur, double FL) ;
	double FE ;
	int nstep,i ;

	if(argc < 3) {
		fprintf(stderr,"usage: inf_const Fhp a [nstep] \n") ;
		exit(0) ;
	}

	sscanf(argv[1],"%lf",&Fhp) ;
	sscanf(argv[2],"%lf",&a) ;

	if(argc > 3)
		sscanf(argv[3],"%d",&nstep) ;
	else
		nstep = 10 ;

	/* initial guesses, good for Fhp = 4, a = 0 */
	x[0] = 3.5 ;	/* fast point Boyer-Lindquist radius */
	x[1] = -0.33 ;	/* u^r at fast point */
	x[2] = -1.9 ;	/* angular momentum flux FL */

	/* initial values for Fhp  and a for which the guesses
	   are good */
	Fhpi = 4. ;
	ai = 0. ;

	/* endpoint of path in Fhp-a plane */
	Fhpf = Fhp ;
	af = a ;

	dFhp = Fhpf - Fhpi ;
	da = af - ai ;

	
	for(i = 0 ; i <= 3*nstep ; i++) {
		if(i < nstep)		/* step in a to 0.5 */
			a = 0.5*i/nstep ;
		else if(i < 2*nstep)	/* step in Fhp */
			a = 0.5 ;
		else			/* step rest of the way in a */
			a = 0.5 + (af - 0.5)*(i - 2*nstep)/nstep ;

		if(i < nstep)		/* step in a to 0.5 */
			Fhp = 4. ;
		else if(i < 2*nstep)	/* step in Fhp */
				Fhp = 4. + (Fhpf - 4.)*(i - nstep)/nstep ;
		else			/* step rest of the way in a */
			Fhp = Fhpf ;

		/* calculate rotation frequency of marginally
		   stable orbit; store in W */
		wmsocalc() ;

		/* solve for fast point, fast point u^r, FL such
		   that FE(fast point) = FE(marginally stable orbit) */
		newt(x-1, 3, &check, vecfunc) ;
	}

	fprintf(stderr,"r_{mso}:     %g \n",rmso) ;
	fprintf(stderr,"r_{horizon}: %g \n",1. + sqrt(1. - a*a)) ;
	fprintf(stderr,"W_{mso} =    %g \n",W) ;

	fprintf(stdout,"fast point radius: %g\n",x[0]) ;
	fprintf(stdout,"u^r at fast point: %g\n\n",x[1]) ;

	FE = efl(rmso,0.,x[2]) ;
	fprintf(stdout,"Use last line as argument to inf_solve.\n") ;
	fprintf(stdout,"Fhp a FL FE:\n") ;
	fprintf(stdout,"%12.8g %12.8g %12.8g %12.8g\n",Fhp,a,x[2],FE) ;

	return(0) ;
}

void vecfunc(int n, double *x, double *fvec)
{
	double rf,uf,FL,Ef,Es ;
	double dedr(double r, double ur, double FL) ;
	double dedur(double r, double ur, double FL) ;
	double efl(double r, double ur, double FL) ;
	
	rf = x[1] ;
	uf = x[2] ;
	FL = x[3] ;

	/* eqtn. 2: fast velocity condition */
	fvec[1] = dedur(rf,uf,FL) ;

	/* eqtn. 4: dE/dr = 0 at fast point */
	fvec[2] = dedr(rf,uf,FL) ;

	/* eqtn. 5: E(fast) = E(slow) */
	Es = efl(rmso,0.,FL) ;
	Ef = efl(rf,uf,FL) ;
	fvec[3] = Ef - Es ;
}

/* calculate \partial (energy flux)/\partial r */
double dedr(double r, double ur, double FL) 
{
	double del ;
	double efl(double r, double ur, double FL) ;

	del = (efl(r*(1. + DEL), ur, FL) - efl(r*(1. - DEL), ur, FL))/
		(efl(r, ur, FL)*2.*DEL) ;

	return(del) ;
}

/* calculate \partial (energy flux)/\partial u^r */
double dedur(double r, double ur, double FL) 
{
	double del ;
	double efl(double r, double ur, double FL) ;

	del = (efl(r, ur*(1. + DEL), FL) - efl(r, ur*(1. - DEL), FL))/
		(efl(r, ur, FL)*2.*DEL) ;

	return(del) ;
}

/* calculate radial energy flux T^r_t as a function of
   radius, u^r, and angular momentum flux T^r_\phi */
double efl(double r,double ur, double FL)
{
	double rho,Fth ;
	double c1,c2,c3,c4,c5,c6,c7 ;
	double d1,d2,d3 ;
	double DD,AA ;
	double up,ut,FE ;
	double Mar,disc,sgn ;

	/* special case for zero velocity */
	if(fabs(ur) < 1.e-12) {
		FE = FM*sqrt(1. - 3./r + 2.*a/pow(r,1.5))/(1 + a/pow(r,1.5))
			+ FL/(pow(r,1.5) + a) ;

		return(FE) ;
	}

	rho = FM/(2.*M_PI*ur*r*r) ;

	Fth = W*Fhp ;

	/* some auxiliary coefficients */
	AA = 1. + a*a/(r*r) + 2.*a*a/(r*r*r) ;
	DD = 1. - 2./r + a*a/(r*r) ;
	c1 = FL ;
	c2 = DD*Fhp*Fhp/(2.*ur) - 2.*AA*rho*M_PI*r*r*r*r*ur ;
	c3 = 4.*a*rho*M_PI*r*ur - DD*Fhp*Fth/(2.*ur) ;
	c4 = 1. + ur*ur/DD ;
	c5 = AA*r*r ;
	c6 = -4.*a/r ;
	c7 = (-1. + 2./r) ;

	d1 = c4 + c1*c1*c5/(c2*c2) ;
	d2 = 2.*c1*c3*c5/(c2*c2) - c1*c6/c2 ;
	d3 = c3*c3*c5/(c2*c2) - c3*c6/c2 + c7 ;

	/* Alfven Mach number */
	Mar = ur*ur*4.*AA*rho*M_PI*pow(r,4)/(DD*Fhp*Fhp) ;

	disc = d2*d2 - 4.*d1*d3 ;

	/* find ut; sgn is the sign of solution, which should
	   have one value inside the Alfven point and another
	   outside */
	if(Mar > 1.) sgn = -1. ;
	else sgn = 1. ;

	if(disc > 0.)
		ut = (-2.*c1*c3*c5 + c1*c2*c6 + sgn*c2*sqrt(
			-4.*c3*c3*c4*c5 + 4.*c2*c3*c4*c6 - 4.*c2*c2*c4*c7 +
			c1*c1*(c6*c6 - 4.*c5*c7)))/(2.*(c3*c3*c5 - 
			c2*c3*c6 + c2*c2*c7)) ;
	else {
		fprintf(stderr,"error in efl\n") ;
		exit(1) ;
	}


	up = -(c1 + c3*ut)/c2 ;

	FE = 2.*M_PI*r*r*(
		DD*Fth*(-Fhp*up + Fth*ut)/(4.*M_PI*r*r*ur)
		+ rho*ur*(2.*a*up/r + (1. - 2./r)*ut)) ;

	return(FE) ;
}

/* calculate orbital frequency of marginally stable orbit */

void wmsocalc()
{
	double Z1,Z2 ;

	/* find radius of mso */
	Z1 = 1. + pow(1. - a*a,1./3.)*(pow(1. + a,1./3.) +
		pow(1. - a, 1./3.)) ;
	Z2 = sqrt(3.*a*a + Z1*Z1) ;
	rmso = 3. + Z2 - sqrt((3. - Z1)*(3. + Z1 + 2.*Z2)) ;

	/* now calculate W for prograde orbit */
	W = 1./(pow(rmso,1.5) + a) ;
}
