#ifndef DTLZ_H
#define DTLZ_H

namespace DTLZ{

	void DTLZ_Vars_Allocate(double *lowBound,double *uppBound,const unsigned int nvar);

	void DTLZ1(double *x, double *f, const unsigned int nvar,const unsigned int nobj);
	void DTLZ2(double *x, double *f, const unsigned int nvar,const unsigned int nobj);
	void DTLZ3(double *x, double *f, const unsigned int nvar,const unsigned int nobj);
	void DTLZ4(double *x, double *f, const unsigned int nvar,const unsigned int nobj);
	void DTLZ5(double *x, double *f, const unsigned int nvar,const unsigned int nobj);
	void DTLZ6(double *x, double *f, const unsigned int nvar,const unsigned int nobj);
	void DTLZ7(double *x, double *f, const unsigned int nvar,const unsigned int nobj);
	//Convex
	void Convex_DTLZ2(double *x, double *f, const unsigned int nvar,const unsigned int nobj);

	//inverted
	void I_DTLZ1(double *x, double *f, const unsigned int nvar, const unsigned int nobj);

	//constrained
	void C1_DTLZ1(double *x, double *f,double *c,const unsigned int nvar,const unsigned int nobj);
	void C1_DTLZ3(double *x, double *f,double *c,const unsigned int nvar,const unsigned int nobj);

	void C2_DTLZ2(double *x, double *f,double *c,const unsigned int nvar,const unsigned int nobj);
	void C2_Convex_DTLZ2(double *x, double *f,double *c,const unsigned int nvar,const unsigned int nobj);

	void C3_DTLZ1(double *x, double *f,double *c,const unsigned int nvar,const unsigned int nobj);
	void C3_DTLZ4(double *x, double *f,double *c,const unsigned int nvar,const unsigned int nobj);

};

#endif