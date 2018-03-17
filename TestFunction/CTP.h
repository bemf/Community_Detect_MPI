#ifndef CTP_H
#define CTP_H

namespace CTP{

	void CTP_Vars_Allocate(double * lowBound,double *uppBound);

	double g_ctp(double *x);
	void calc_cv_ctp1(double *constain,double *f, double g,double a1,double a2,double b1,double b2);
	void calc_cv_ctp(double *constain,double *f,double g,double theta,double a,double b,double c,double d,double e);

	void CTP1(double *x, double *f,double *c);
	void CTP2(double *x, double *f,double *c);
	void CTP3(double *x, double *f,double *c);
	void CTP4(double *x, double *f,double *c);
	void CTP5(double *x, double *f,double *c);
	void CTP6(double *x, double *f,double *c);
	void CTP7(double *x, double *f,double *c);
};

#endif