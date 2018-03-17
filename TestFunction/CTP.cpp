#include "CTP.h"
#include <math.h>

namespace CTP{
	#define PI  3.1415926535897932384626433832795
	
	void CTP_Vars_Allocate(double * lowBound,double *uppBound)
	{
		lowBound[0] = 0.0;
		lowBound[1] = 0.0;
		uppBound[0] = 1.0;
		uppBound[1] = 5.0;
	}

	double g_ctp(double *x)
	{
		return 1.0 + x[1];
	}
	void  calc_cv_ctp1(double *constain,double *f, double g,double a1,double a2,double b1,double b2)
	{
		constain[0]= f[1] - a1 * exp(-b1 * f[0]);
		constain[1] = f[1] - a2 * exp(-b2 * f[0]);
	}
	void calc_cv_ctp(double *constain,double *f,double g,double theta,double a,double b,double c,double d,double e)
	{
		constain[0] = cos(theta)*(f[1]-e) - sin(theta)*f[0] - a*pow(fabs(sin(b*PI*pow(sin(theta)*(f[1] - e)+cos(theta)*f[0],c))),d);
	}

	void CTP1(double *x, double *f,double *c)
	{
		f[0] = x[0];
        double g = g_ctp(x);
        f[1] = g * exp(-f[0]/g);
        calc_cv_ctp1(c,f,g,0.858,0.728,0.541,0.295);
	}

	void CTP2(double *x, double *f,double *c)
	{
		f[0] = x[0];
        double g = g_ctp(x);
        f[1] = g*(1.0 - sqrt(f[0]/g));
        calc_cv_ctp(c,f,g,-0.2*PI,0.2,10,1,6,1);
	}

	void CTP3(double *x, double *f,double *c)
	{
		 f[0] = x[0];
        double g = g_ctp(x);
        f[1] = g*(1.0 - sqrt(f[0]/g));
        calc_cv_ctp(c,f,g,-0.2*PI,0.1,10,1,0.5,1);
	}

	void CTP4(double *x, double *f,double *c)
	{
		 f[0] = x[0];
        double g = g_ctp(x);
        f[1] = g*(1.0 - sqrt(f[0]/g));
        calc_cv_ctp(c,f,g,-0.2*PI,0.75,10,1,0.5,1);
	}

	void CTP5(double *x, double *f,double *c)
	{
		 f[0] = x[0];
        double g = g_ctp(x);
        f[1] = g*(1.0 - sqrt(f[0]/g));
        calc_cv_ctp(c,f,g,-0.2*PI,0.75,10,2,0.5,1);
	}

	void CTP6(double *x, double *f,double *c)
	{
		 f[0] = x[0];
        double g = g_ctp(x);
        f[1] = g*(1.0 - sqrt(f[0]/g));
        calc_cv_ctp(c,f,g,0.1*PI,40,0.5,1,2,-2);
	}

	void CTP7(double *x, double *f,double *c)
	{
		 f[0] = x[0];
        double g = g_ctp(x);
        f[1] = g*(1.0 - sqrt(f[0]/g));
        calc_cv_ctp(c,f,g,-0.05*PI,40,5,1,6,0);
	}
}