// Functions to sample from Poisson, binomial
// and normal distributions (from Numerical Recipes in C)

#include "functions.h"
#include <cmath>
using namespace std;

extern MTRand rnd;


double gammln(const double xx)
{
	int j;
	double x,y,tmp,ser;
	static const double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,0.1208650973866179e-2,
		-0.5395239384953e-5};
    
	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<6;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

double poisdev(const double xm)
{
	const double PI=3.141592653589793238;
	static double sq,alxm,g,oldm=(-1.0);
	double em,t,y;
    
	if (xm < 12.0) {
		if (xm != oldm) {
			oldm=xm;
			g=exp(-xm);
		}
		em = -1;
		t=1.0;
		do {
			++em;
			t *= rnd.rand();;
		} while (t > g);
	} else {
		if (xm != oldm) {
			oldm=xm;
			sq=sqrt(2.0*xm);
			alxm=log(xm);
			g=xm*alxm-gammln(xm+1.0);
		}
		do {
			do {
				y=tan(PI*rnd.rand());
				em=sq*y+xm;
			} while (em < 0.0);
            em=floor(em);
            t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
		} while (rnd.rand() > t);
	}
	return em;
}


double binldev(const double pp, const int n)
{
	const double PI=3.141592653589793238;
	int j;
	static int nold=(-1);
	double am,em,g,angle,p,bnl,sq,t,y;
	static double pold=(-1.0),pc,plog,pclog,en,oldg;
    
	p=(pp <= 0.5 ? pp : 1.0-pp);
	am=n*p;
	if (n < 25) {
		bnl=0.0;
		for (j=0;j<n;j++)
			if (rnd.rand() < p) ++bnl;
	} else if (am < 1.0) {
		g=exp(-am);
		t=1.0;
		for (j=0;j<=n;j++) {
			t *= rnd.rand();
			if (t < g) break;
		}
		bnl=(j <= n ? j : n);
	} else {
		if (n != nold) {
			en=n;
			oldg=gammln(en+1.0);
			nold=n;
		} if (p != pold) {
			pc=1.0-p;
			plog=log(p);
			pclog=log(pc);
			pold=p;
		}
		sq=sqrt(2.0*am*pc);
		do {
			do {
				angle=PI*rnd.rand();
				y=tan(angle);
				em=sq*y+am;
			} while (em < 0.0 || em >= (en+1.0));
			em=floor(em);
			t=1.2*sq*(1.0+y*y)*exp(oldg-gammln(em+1.0)
                                   -gammln(en-em+1.0)+em*plog+(en-em)*pclog);
		} while (rnd.rand() > t);
		bnl=em;
	}
	if (p != pp) bnl=n-bnl;
	return bnl;
}


// draw from a normal distribution:

double gasdev()
{
	static int iset=0;
	static double gset;
	double fac,rsq,v1,v2;
	
	//	if (idum < 0) iset=0;
	if (iset == 0) {
		do {
			v1=2.0*rnd.rand()-1.0;
			v2=2.0*rnd.rand()-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
};
