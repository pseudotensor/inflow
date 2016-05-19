

/* 

inf_solv: find inflow solution.  Use this after using
inf_const to find appropriate constants to use as
arguments.

arguments:

Fhp: F_{\theta\phi} component of maxwell tensor, basically
	the radial magnetic field strength 
 
a: spin parameter for black hole 

FL: radial angular momentum flux

FE: radial energy flux 

N: number of points in output

Output lines are
r, rho, u^r, u^\phi, B^r, B^\phi, Alfven Mach no.
Here B^i is F^{*it}.

Output can be adjusted to give constrained transport
variables F_{ij}, which are related to the B^i's by
a factor of \sqrt{-g}.

NOTES

*solution can be ratty for high values of 
magnetic field Fhp, or for very high resolution,
where the solver spends a lot of time near the
Alfven point. Values that are known to work
are Fhp = 0.5, a = 0.5.

CFG 3-20-99
updated 7-19-02 for testing grmhd code.

*/

#include <stdio.h>
#include <math.h>

/** global variables **/

double FM = -1. ;	/* mass flux */
double FL ;		/* angular momentum flux */
double FE ;		/* energy flux */
double r ;		/* current radius */
double GM = 1. ;	/* duh */
double Fhp ;		/* theta-phi component of Maxwell; reduces to -r^2 Br
			   in nonrelativistic limit */
double W ;		/* "orbital frequency" of magnetic field, chosen
			   to be orbital frequency at marginally stable
			   orbit */
double a ;		/* black hole spin parameter, |a| < 1.
			   currently prograde orbits around prograde hole
			   are assumed. */
double rmso ;		/* radius of marginally stable orbit */
double rh ;		/* horizon radius */
double disc ;
double sgn ;

double DTOL = 1.e-4 ;
double DEL = 2.e-5 ;	/* optimal choice for centered difference 
			   numerical derivatives */

main(int argc, char *argv[])
{
	double efl(double ur) ;
	void varcalc(double ur,double *up,double *ut,double *rho,
		double *l,double *E,double *Fth,double *Frh,
		double *Mar,double *Maf,double *ed)  ;
	double find_root(double (*func)(), double guess, double tol) ;
	double ur,uh,up,ut,rho,l,E,Fth,Frh,Mar,Maf,ed ;
	double rh,dr,urguess,urlast,urlastlast,rs,us,r0 ;
	double g,Br,Bh,Bp,rrat,srat ;
	double DD,Trpfl,Trpem ;
	double lr,lr0,lr1,dlr,rp ;
	int i,nstep,N ;
	void wmsocalc() ;

	if(argc < 5) {
		fprintf(stderr,"usage: inf_solv Fhp a FL FE N \n") ;
		exit(0) ;
	}

	sscanf(argv[1],"%lf",&Fhp) ;
	sscanf(argv[2],"%lf",&a) ;
	sscanf(argv[3],"%lf",&FL) ;
	sscanf(argv[4],"%lf",&FE) ;
	sscanf(argv[5],"%d",&N) ;

	/* find rotation frequency of marginally stable orbit */
        wmsocalc() ;

	/* find horizon radius */
	rh = 1. + sqrt(1. - a*a) ;

	sgn = -1. ;	/* don't fiddle w/ this */

	r0 = 1.02*rh ;
	r = r0 ;
	urguess = find_root(efl,-0.4,1.e-9) ;
	urlast = urguess ;
	
	lr0 = log(r0) ;
	lr1 = log(0.98*rmso) ;
	dlr = (lr1 - lr0)/(N+2.) ;

	/* use a logarithmic radial grid */
	for(lr = lr0 + 0.5*dlr ; lr < lr1 ; lr += dlr ) {
		r = exp(lr) ;
		rp = r0*exp((i+0.5)*dlr) ;

		/* use this to set sign in ut solution */
		if(urguess != urlast) {
			varcalc(urguess,&up,&ut,&rho,&l,&E,&Fth,&Frh,&Mar,&Maf,&ed) ;
			if(Mar < 1.) sgn = 1. ;
			else sgn = -1. ;
		}
		ur = find_root(efl,urguess,1.e-10) ;

		/* given u^r, find all other quantities associated with
		   solution */
		varcalc(ur,
			&up,	/* u^\phi */
			&ut,	/* u^t */
			&rho,	/* density */
			&l,	/* u_\phi */
			&E,	/* u_t (or -u_t...) */
			&Fth,	/* F_{t\theta} */
			&Frh,	/* F_{r\theta} */
			&Mar,	/* fast mach no. */
			&Maf,	/* Alfven mach no. */
			&ed	/* electromagnetic energy density */
			) ;

		/* output results */
		uh = 0. ;
		fprintf(stdout,"%15.10g %15.10g %15.10g %15.10g ",
				r,rho,ur,up) ;

		/* calculate field components */
		g = r*r ;
		Br = Fhp/g/sqrt(4.*M_PI) ;
		Bp = Frh/g/sqrt(4.*M_PI) ;
		fprintf(stdout,"%15.10g %15.10g\n",Br,Bp) ;

		/* Hawley-deVilliers field output */
		/* CAUTION: factors of sqrt(4 Pi) may be present--
		   these components have not had the sqrt(4 Pi)
		   absorbed into them yet */
		/*
		fprintf(stdout,"%15.10g %15.10g ",Fhp,Frh) ;
		*/

		/* other quantities of interest */
		DD = 1. - 2./r + a*a/(r*r) ;
		Trpfl = rho*ur*l ;
		Trpem = DD*Fhp*Frh/(4.*M_PI*r*r) ;
		/*
		fprintf(stdout,"%15.10g %15.10g %15.10g\n",Mar,Trpfl,Trpem) ;
		*/

		/* guess at next root, to keep on correct branch */
		urguess = ur + (ur - urlast) ;
		urlast = ur ;
	}
}

double efl(double ur)
{
	double rho,Fth ;
	double c1,c2,c3,c4,c5,c6,c7 ;
	double d1,d2,d3 ;
	double DD,AA ;
	double up,ut,FEp ;
	double Mar ;

	rho = FM/(2.*M_PI*ur*r*r) ;

	/* consequence of tying inflow to disk */
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

	disc = d2*d2 - 4.*d1*d3 ;

	/* find ut; sgn chooses solution.  Should 
	   have one sign inside the Alfven point and another
	   outside */
	if(disc > 0.) {
		ut = (-2.*c1*c3*c5 + c1*c2*c6 + sgn*c2*sqrt(
			-4.*c3*c3*c4*c5 + 4.*c2*c3*c4*c6 - 4.*c2*c2*c4*c7 +
			c1*c1*(c6*c6 - 4.*c5*c7)))/(2.*(c3*c3*c5 - c2*c3*c6 + c2*c2*c7)) ;
	}
	else {
                ut = (-2.*c1*c3*c5 + c1*c2*c6 - sgn*c2*sqrt(-(
                        -4.*c3*c3*c4*c5 + 4.*c2*c3*c4*c6 - 4.*c2*c2*c4*c7 +
                        c1*c1*(c6*c6 - 4.*c5*c7))))/(2.*(c3*c3*c5 - c2*c3*c6 + c2*c2*c7)) ;
	}

	up = -(c1 + c3*ut)/c2 ;

	FEp = 2.*M_PI*r*r*(
		DD*Fth*(-Fhp*up + Fth*ut)/(4.*M_PI*r*r*ur)
		+ rho*ur*(2.*a*up/r + (1. - 2./r)*ut)) ;

	return(FEp - FE) ;
}

void varcalc(double ur,double *up,double *ut,double *rho,
	double *l,double *E,double *Fth,double *Frh,
	double *Mar,double *Maf, double *ed)
{
	double Up,Br,rhow,k0,Uasq,Bp,f1,Ufsq ;
	double AA,DD,c1,c2,c3,c4,c5,c6,c7,d1,d2,d3,disc ;
	double ed_calc(double a, double Fhp, double r, double ur,
		double up, double ut, double Fth, double Frh) ;

	*rho = FM/(2.*M_PI*ur*r*r) ;

	*Fth = W*Fhp ;

	/* some auxiliary coefficients */
	AA = 1. + a*a/(r*r) + 2.*a*a/(r*r*r) ;
	DD = 1. - 2./r + a*a/(r*r) ;
	c1 = FL ;
	c2 = DD*Fhp*Fhp/(2.*ur) - 2.*AA*(*rho)*M_PI*r*r*r*r*ur ;
	c3 = 4.*a*(*rho)*M_PI*r*ur - DD*Fhp*(*Fth)/(2.*ur) ;
	c4 = 1. + ur*ur/DD ;
	c5 = AA*r*r ;
	c6 = -4.*a/r ;
	c7 = (-1. + 2./r) ;

	d1 = c4 + c1*c1*c5/(c2*c2) ;
	d2 = 2.*c1*c3*c5/(c2*c2) - c1*c6/c2 ;
	d3 = c3*c3*c5/(c2*c2) - c3*c6/c2 + c7 ;

	disc = d2*d2 - 4.*d1*d3 ;

	if(disc > 0.)
		*ut = (-2.*c1*c3*c5 + c1*c2*c6 + sgn*c2*sqrt(
			-4.*c3*c3*c4*c5 + 4.*c2*c3*c4*c6 - 4.*c2*c2*c4*c7 +
			c1*c1*(c6*c6 - 4.*c5*c7)))/(2.*(c3*c3*c5 - c2*c3*c6 + c2*c2*c7)) ;
	else
		*ut = (-2.*c1*c3*c5 + c1*c2*c6 - sgn*c2*sqrt(-(
			-4.*c3*c3*c4*c5 + 4.*c2*c3*c4*c6 - 4.*c2*c2*c4*c7 +
			c1*c1*(c6*c6 - 4.*c5*c7))))/(2.*(c3*c3*c5 - c2*c3*c6 + c2*c2*c7)) ;

	*up = -(c1 + c3*(*ut))/c2 ;

	*l = AA*r*r*(*up) - 2.*a*(*ut)/r ;

	*E = (1. - 2./r)*(*ut) + 2.*a*(*up)/r ;

	*Frh = (Fhp*(*up) - (*Fth)*(*ut))/ur ;

	/* following Takahashi et al.'s appendix */
	Up = ur/sqrt(DD) ;
	Br = Fhp/(r*r*sqrt(DD)) ;
	rhow = r*sqrt(DD) ;
	k0 = 1. - 2./r + 4.*a*W/r - AA*r*r*W*W ;
	Uasq = ( Br*Br*k0 / ( 4.*M_PI*(*rho) ) ) ;
	*Mar = sqrt(Up*Up/Uasq) ;
	Bp = -DD*(*Frh) ;

	/* slow and fast speed */
	Ufsq = Uasq + Bp*Bp/(4.*M_PI*(*rho)*rhow*rhow) ;

	*Maf = sqrt(Up*Up/Ufsq) ;

	*ed = ed_calc(a,Fhp,r,ur,*up,*ut,*Fth,*Frh) ;
}

#define GENXTOL	1.e-10

/* silly little root finder; given all the constants of the
 * solution, find u^r */
double find_root(double (*func)(), double guess, double tol)
{
        double x1, x2, del, y1, y2, xnew, ynew, max ;

        x1 = guess ;
        y1 = (*func)(x1) ;
        del = -1.e-5*guess ;
        x2 = x1 + del ;
        y2 = (*func)(x2) ;

        xnew = (x2*y1 - x1*y2)/(y1 - y2) ;
        ynew = (*func)(xnew) ;

        while(fabs(ynew) > tol && fabs(x1 - x2) > GENXTOL*fabs(x1)) {
                /* revert to stepping if estimate is too far away */
                max = 2.*fabs(x1 - x2) ;
                if(fabs(xnew - x1) > max && fabs(xnew - x2) > max) {
                        /* x2 is closer */
                        if(fabs(xnew - x1) > fabs(xnew - x2)) {
                                xnew = x2 + max*(xnew - x2)/fabs(xnew - x2) ;
                                x1 = xnew ;
                        }
                        /* x1 is closer */
                        else {
                                xnew = x1 + max*(xnew - x1)/fabs(xnew - x1) ;
                                x2 = x1 ;
                                x1 = xnew ;
                        }

                }
                else {
                        /* first check to see if root is bracketed */
                        if(y2*ynew < 0.) {
                                x1 = xnew ;
                        }
                        else if(y1*ynew < 0.) {
                                x2 = x1 ;
                                x1 = xnew ;
                        }       
                        else {
                                /* x2 is closer */ 
                                if(fabs(xnew - x1) > fabs(xnew - x2)) {
                                        xnew = x2 + max*(xnew - x2)/fabs(xnew - x2) ;
                                        x1 = xnew ;
                                }
                                /* x1 is closer */
                                else {
                                        xnew = x1 + max*(xnew - x1)/fabs(xnew - x1) ;
                                        x2 = x1 ; 
                                        x1 = xnew ;
                                }
                        }
                }


                y1 = (*func)(x1) ;
                y2 = (*func)(x2) ;

                xnew = (x2*y1 - x1*y2)/(y1 - y2) ;
                ynew = (*func)(xnew) ;

        }

        return(xnew) ;
}

void wmsocalc()
{
        double Z1,Z2 ;

        /* first find rmso */
        Z1 = 1. + pow(1. - a*a,1./3.)*(pow(1. + a,1./3.) +
                pow(1. - a, 1./3.)) ;
        Z2 = sqrt(3.*a*a + Z1*Z1) ;
        rmso = 3. + Z2 - sqrt((3. - Z1)*(3. + Z1 + 2.*Z2)) ;
        fprintf(stderr,"rmso: %g \n",rmso) ;

        /* now calculate W for prograde orbit */
        W = 1./(pow(rmso,1.5) + a) ;
        fprintf(stderr,"W = %g\n",W) ;
}

/* calculate electromagnetic energy density in comoving frame */
double ed_calc(double a, double Fhp, double r, double ur,
	double up, double ut, double Fth, double Frh) 
{
	double r2,r3,r4,r5,d,DD,AA,Frh2,ed ;

	r2 = r*r ;
	r3 = r2*r ;
	r4 = r3*r ;
	r5 = r4*r ;

	d = -1./(2.*M_PI*r2*ur) ;

	DD = 1. - 2./r + pow(a/r,2) ;
	AA = 1. + pow(a/r,2) + 2.*pow(a,2)/pow(r,3) ;
	Frh2 = Frh*Frh ;
	
	/* there are smarter ways to do this */
	ed = ((AA*AA*r4*(-((-(Fhp*
	                          ((-2.*a*Fth)/(DD*r5) - (Fhp*(1. - 2./r))/(DD*r4))) 
	                        + Fhp*((2.*a*Fth)/(DD*r5) + 
	                          (Fhp*(1. - 2./r))/(DD*r4)) + 
	                       Fth*((2.*a*Fhp)/(DD*r5) - (AA*Fth)/(DD*r2)) - 
	                       Fth*((-2.*a*Fhp)/(DD*r5) + (AA*Fth)/(DD*r2)) + 
	                       (2.*DD*Frh2)/r2)*(1. - 2./r))/(4.*DD*r2) + 
	                 pow((-2.*a*Fth)/(DD*r5) - (Fhp*(1. - 2./r))/(DD*r4),2)*r2))/
	             (4.*M_PI) - (a*AA*r*((a*
	                    (-(Fhp*((-2.*a*Fth)/(DD*r5) - 
	                           (Fhp*(1. - 2./r))/(DD*r4))) + 
	                      Fhp*((2.*a*Fth)/(DD*r5) + (Fhp*(1. - 2./r))/(DD*r4)) + 
	                      Fth*((2.*a*Fhp)/(DD*r5) - (AA*Fth)/(DD*r2)) - 
	                      Fth*((-2.*a*Fhp)/(DD*r5) + (AA*Fth)/(DD*r2)) + 
	                      (2.*DD*Frh2)/r2))/(2.*DD*r3) + 
	                 ((-2.*a*Fth)/(DD*r5) - (Fhp*(1. - 2./r))/(DD*r4))*
	                  ((2.*a*Fhp)/(DD*r5) - (AA*Fth)/(DD*r2))*r2))/M_PI + 
	            (a*a*((AA*(-(Fhp*((-2.*a*Fth)/(DD*r5) - 
	                           (Fhp*(1. - 2./r))/(DD*r4))) + 
	                      Fhp*((2.*a*Fth)/(DD*r5) + (Fhp*(1. - 2./r))/(DD*r4)) + 
	                      Fth*((2.*a*Fhp)/(DD*r5) - (AA*Fth)/(DD*r2)) - 
	                      Fth*((-2.*a*Fhp)/(DD*r5) + (AA*Fth)/(DD*r2)) + 
	                      (2.*DD*Frh2)/r2))/(4.*DD) + 
	                 pow((2.*a*Fhp)/(DD*r5) - (AA*Fth)/(DD*r2),2)*r2))/
	             (M_PI*r2))*up*up + 
	         2.*(-(a*Frh*((2.*a*Fhp)/(DD*r5) - (AA*Fth)/(DD*r2)))/(2.*M_PI*r) + 
	            (AA*Frh*((-2.*a*Fth)/(DD*r5) - (Fhp*(1. - 2./r))/(DD*r4))*r2)/
	             (4.*M_PI))*up*ur + ((-(DD*
	                  (-(Fhp*((-2.*a*Fth)/(DD*r5) - 
	                         (Fhp*(1. - 2./r))/(DD*r4))) + 
	                    Fhp*((2.*a*Fth)/(DD*r5) + (Fhp*(1. - 2./r))/(DD*r4)) + 
	                    Fth*((2.*a*Fhp)/(DD*r5) - (AA*Fth)/(DD*r2)) - 
	                    Fth*((-2.*a*Fhp)/(DD*r5) + (AA*Fth)/(DD*r2)) + 
	                    (2.*DD*Frh2)/r2))/4. + (DD*DD*Frh2)/r2)*ur*ur)/
	          (4.*DD*DD*M_PI) + 2.*(-(a*AA*r*
	                (-((-(Fhp*((-2.*a*Fth)/(DD*r5) - 
	                             (Fhp*(1. - 2./r))/(DD*r4))) + 
	                        Fhp*((2.*a*Fth)/(DD*r5) + 
	                           (Fhp*(1. - 2./r))/(DD*r4)) + 
	                        Fth*((2.*a*Fhp)/(DD*r5) - (AA*Fth)/(DD*r2)) - 
	                        Fth*((-2.*a*Fhp)/(DD*r5) + (AA*Fth)/(DD*r2)) + 
	                        (2.*DD*Frh2)/r2)*(1. - 2./r))/(4.*DD*r2) + 
	                  pow((-2.*a*Fth)/(DD*r5) - (Fhp*(1. - 2./r))/(DD*r4),2)*r2))
	              /(2.*M_PI) + (a*a*((a*
	                    (-(Fhp*((-2.*a*Fth)/(DD*r5) - 
	                           (Fhp*(1. - 2./r))/(DD*r4))) + 
	                      Fhp*((2.*a*Fth)/(DD*r5) + (Fhp*(1. - 2./r))/(DD*r4)) + 
	                      Fth*((2.*a*Fhp)/(DD*r5) - (AA*Fth)/(DD*r2)) - 
	                      Fth*((-2.*a*Fhp)/(DD*r5) + (AA*Fth)/(DD*r2)) + 
	                      (2.*DD*Frh2)/r2))/(2.*DD*r3) + 
	                 ((-2.*a*Fth)/(DD*r5) - (Fhp*(1. - 2./r))/(DD*r4))*
	                  ((2.*a*Fhp)/(DD*r5) - (AA*Fth)/(DD*r2))*r2))/(M_PI*r2)
	             + (AA*(-1. + 2./r)*r2*
	               ((a*(-(Fhp*((-2.*a*Fth)/(DD*r5) - 
	                           (Fhp*(1. - 2./r))/(DD*r4))) + 
	                      Fhp*((2.*a*Fth)/(DD*r5) + (Fhp*(1. - 2./r))/(DD*r4)) + 
	                      Fth*((2.*a*Fhp)/(DD*r5) - (AA*Fth)/(DD*r2)) - 
	                      Fth*((-2.*a*Fhp)/(DD*r5) + (AA*Fth)/(DD*r2)) + 
	                      (2.*DD*Frh2)/r2))/(2.*DD*r3) + 
	                 ((-2.*a*Fth)/(DD*r5) - (Fhp*(1. - 2./r))/(DD*r4))* 
	                  ((2.*a*Fhp)/(DD*r5) - (AA*Fth)/(DD*r2))*r2))/(4.*M_PI) - 
	            (a*(-1. + 2./r)*((AA*(-(Fhp*
	                         ((-2.*a*Fth)/(DD*r5) - (Fhp*(1. - 2./r))/(DD*r4))) 
	                       + Fhp*((2.*a*Fth)/(DD*r5) + 
	                         (Fhp*(1. - 2./r))/(DD*r4)) + 
	                      Fth*((2.*a*Fhp)/(DD*r5) - (AA*Fth)/(DD*r2)) - 
	                      Fth*((-2.*a*Fhp)/(DD*r5) + (AA*Fth)/(DD*r2)) + 
	                      (2.*DD*Frh2)/r2))/(4.*DD) + 
	                 pow((2.*a*Fhp)/(DD*r5) - (AA*Fth)/(DD*r2),2)*r2))/
	             (2.*M_PI*r))*up*ut + 
	         2.*((Frh*((2.*a*Fhp)/(DD*r5) - (AA*Fth)/(DD*r2))*(-1. + 2./r))/
	             (4.*M_PI) - (a*Frh*((-2.*a*Fth)/(DD*r5) - 
	                 (Fhp*(1. - 2./r))/(DD*r4)))/(2.*M_PI*r))*ur*ut + 
	         ((a*a*(-((-(Fhp*((-2.*a*Fth)/(DD*r5) - 
	                            (Fhp*(1. - 2./r))/(DD*r4))) + 
	                       Fhp*((2.*a*Fth)/(DD*r5) + 
	                          (Fhp*(1. - 2./r))/(DD*r4)) + 
	                       Fth*((2.*a*Fhp)/(DD*r5) - (AA*Fth)/(DD*r2)) - 
	                       Fth*((-2.*a*Fhp)/(DD*r5) + (AA*Fth)/(DD*r2)) + 
	                       (2.*DD*Frh2)/r2)*(1. - 2./r))/(4.*DD*r2) + 
	                 pow((-2.*a*Fth)/(DD*r5) - (Fhp*(1. - 2./r))/(DD*r4),2)*r2))/
	             (M_PI*r2) - (a*(-1. + 2./r)*
	               ((a*(-(Fhp*((-2.*a*Fth)/(DD*r5) - 
	                           (Fhp*(1. - 2./r))/(DD*r4))) + 
	                      Fhp*((2.*a*Fth)/(DD*r5) + (Fhp*(1. - 2./r))/(DD*r4)) + 
	                      Fth*((2.*a*Fhp)/(DD*r5) - (AA*Fth)/(DD*r2)) - 
	                      Fth*((-2.*a*Fhp)/(DD*r5) + (AA*Fth)/(DD*r2)) + 
	                      (2.*DD*Frh2)/r2))/(2.*DD*r3) + 
	                 ((-2.*a*Fth)/(DD*r5) - (Fhp*(1. - 2./r))/(DD*r4))*
	                  ((2.*a*Fhp)/(DD*r5) - (AA*Fth)/(DD*r2))*r2))/(M_PI*r) + 
	            (pow(-1. + 2./r,2)*((AA*
	                    (-(Fhp*((-2.*a*Fth)/(DD*r5) - 
	                           (Fhp*(1. - 2./r))/(DD*r4))) + 
	                      Fhp*((2.*a*Fth)/(DD*r5) + (Fhp*(1. - 2./r))/(DD*r4)) + 
	                      Fth*((2.*a*Fhp)/(DD*r5) - (AA*Fth)/(DD*r2)) - 
	                      Fth*((-2.*a*Fhp)/(DD*r5) + (AA*Fth)/(DD*r2)) + 
	                      (2.*DD*Frh2)/r2))/(4.*DD) + 
	                 pow((2.*a*Fhp)/(DD*r5) - (AA*Fth)/(DD*r2),2)*r2))/(4.*M_PI))
	           *ut*ut ;

	return(ed) ;
}
