#include "DTLZ.h"
#include <math.h>

namespace DTLZ{
	#define PI  3.1415926535897932384626433832795
	void DTLZ_Vars_Allocate(double *lowBound,double *uppBound,const unsigned int nvar)
	{
		for(int i = 0 ; i < nvar;i++)
		{
			lowBound[i]=0.0;
			uppBound[i]=1.0;
		}
	}

	void DTLZ1(double *x, double *f,const unsigned int nvar,const unsigned int nobj)
	{
		double g = 0;
		int k = nvar - nobj + 1;
		for (int i = nvar - k; i < nvar; i++)
		{
			g += (pow((x[i] - 0.5), 2) - cos(20 * PI * (x[i] - 0.5)));
		}
		g = 100 * (k + g);

		for (int i = 0; i < nobj; i++)
		{
			f[i] = 0.5*(1.0 + g);

			for (int j = 0; j < nobj - i - 1; j++)
			{
				f[i] *= x[j];
			}

			if (i != 0)
			{
				f[i] *= 1 - x[nobj - i - 1];
			}
		}
	}
	void DTLZ2(double *x, double *f, const unsigned int nvar,const unsigned int nobj)
	{
		double g = 0;
		int k = nvar - nobj + 1;

		for (int i = nvar - k; i < nvar; i++)
		{
			g += pow(x[i] - 0.5, 2.0);
		}

		for (int i = 0; i < nobj; i++)
		{
			f[i] = 1.0 + g;

			for (int j = 0; j < nobj - i - 1; j++)
			{
				f[i] *= cos(0.5*PI*x[j]);
			}

			if (i != 0)
			{
				f[i] *= sin(0.5*PI*x[nobj - i - 1]);
			}
		}
	}

	void DTLZ3(double *x, double *f, const unsigned int nvar,const unsigned int nobj)
	{
		double g = 0;
		int k = nvar - nobj + 1;

		for (int i = nvar - k; i < nvar; i++)
		{
			g += pow(x[i] - 0.5, 2.0) - cos(20.0*PI*(x[i] - 0.5));
		}
		g = 100.0 * (k + g);

		for (int i = 0; i < nobj; i++)
		{
			f[i] = 1.0 + g;

			for (int j = 0; j < nobj - i - 1; j++)
			{
				f[i] *= cos(0.5*PI*x[j]);
			}

			if (i != 0)
			{
				f[i] *= sin(0.5*PI*x[nobj - i - 1]);
			}
		}
	}

	void DTLZ4(double *x, double *f, const unsigned int nvar,const unsigned int nobj)
	{
		double alpha = 100.0;
		double g = 0;

		int k = nvar - nobj + 1;

		for (int i = nvar - k; i < nvar; i++)
		{
			g += pow(x[i] - 0.5, 2.0);
		}

		for (int i = 0; i < nobj; i++)
		{
			f[i] = 1.0 + g;

			for (int j = 0; j < nobj - i - 1; j++)
			{
				f[i] *= cos(0.5*PI*pow(x[j], alpha));
			}

			if (i != 0)
			{
				f[i] *= sin(0.5*PI*pow(x[nobj - i - 1], alpha));
			}
		}
	}
	//未检查验证
	void DTLZ5(double *x, double *f, const unsigned int nvar,const unsigned int nobj)
	{
		double g = 0;
		double *theta = new double[nobj];
		int k = nvar - nobj + 1;

		for (int i = nvar - k; i < nvar; i++)
		{
			g += pow(x[i] - 0.5, 2.0);
		}
		double t = PI / (4 * (1 + g));
		theta[0] = 0.5 * x[0] * PI;
		for (int i = 1; i < nobj; i++)
		{
			theta[i] = t * (1 + 2 * g * x[i]);
		}

		for (int i = 0; i < nobj; i++)
		{
			f[i] = 1.0 + g;

			for (int j = 0; j < nobj - i - 1; j++)
			{
				f[i] *= cos(theta[j]);
			}

			if (i != 0)
			{
				f[i] *= sin(theta[nobj - i]);
			}
		}
		delete theta;
	}
	//未检查验证
	void DTLZ6(double *x, double *f, const unsigned int nvar,const unsigned int nobj)
	{
		double g = 0;
		double *theta = new double[nobj];
		int k = nvar - nobj + 1;

		for (int i = nvar - k; i < nvar; i++)
		{
			g += pow(x[i], 0.1);
		}
		double t = PI / (4 * (1 + g));
		theta[0] = 0.5 * x[0] * PI;
		for (int i = 1; i < nobj; i++)
		{
			theta[i] = t * (1 + 2 * g * x[i]);
		}

		for (int i = 0; i < nobj; i++)
		{
			f[i] = 1.0 + g;

			for (int j = 0; j < nobj - i - 1; j++)
			{
				f[i] *= cos(theta[j]);
			}

			if (i != 0)
			{
				f[i] *= sin(theta[nobj - i]);
			}
		}
		delete theta;
	}
	//未检查验证
	void DTLZ7(double *x, double *f, const unsigned int nvar,const unsigned int nobj)
	{
		double g = 0;
		int k = nvar - nobj + 1;

		for (int i = nvar - k; i < nvar; i++)
		{
			g += x[i];
		}
		g = 1.0 + (9.0*g) / k;

		double h = nobj;
		for (int i = 0; i < nobj - 1; i++)
		{
			h -= x[i] / (1.0 + g)*(1.0 + sin(3.0*PI*x[i]));
		}

		for (int i = 0; i < nobj - 1; i++)
		{
			f[i] = x[i];
		}
		f[nobj - 1] = (1.0 + g) * h;
	}

	void Convex_DTLZ2(double *x, double *f, const unsigned int nvar,const unsigned int nobj)
	{
		DTLZ::DTLZ2(x,f,nvar,nobj);
		for (int i = 0; i < nobj - 1; i++)
		{
			f[i] = pow(f[i], 4.0);
		}
		f[nobj - 1] = pow(f[nobj - 1], 2.0);
	}

	void I_DTLZ1(double *x, double *f, const unsigned int nvar, const unsigned int nobj)
	{
		DTLZ::DTLZ1(x, f, nvar, nobj);

		double g = 0;
		int k = nvar - nobj + 1;
		for (int i = nvar - k; i < nvar; i++)
		{
			g += (pow((x[i] - 0.5), 2) - cos(20 * PI * (x[i] - 0.5)));
		}
		g = 100 * (k + g);

		for (int i = 0; i < nobj; i++)
		{
			f[i] = 0.5 * (1 + g) - f[i];
		}

	}

	void C1_DTLZ1(double *x, double *f,double *c,const unsigned int nvar,const unsigned int nobj)
	{
		DTLZ::DTLZ1(x,f,nvar,nobj);
		
		double tmp = 0;
		for(int i =0 ; i< nobj-1;i++)
			tmp += f[i]*2;

		c[0] = 1.0 - f[nobj-1]/0.6 - tmp;
	}

	void C1_DTLZ3(double *x, double *f,double *c,const unsigned int nvar,const unsigned int nobj)
	{
		DTLZ::DTLZ3(x,f,nvar,nobj);
		
		double r[] = { 6, 9, 9, 12.5, 12.5, 12.5, 12.5, 12.5, 15, 15, 15, 15, 15, 15 };
		double tmp1 = 0;
		double tmp2 = 0;
		for(int i = 0 ; i < nobj;i++)
		{
			tmp1 += pow(f[i],2.0);
			tmp2 += pow(f[i], 2.0); 
		}
		c[0] = (tmp1 - 16) * (tmp2 - pow(r[nobj - 2], 2.0));
	}

	void C2_DTLZ2(double *x, double *f, double *c, const unsigned int nvar, const unsigned int nobj)
	{
		DTLZ::DTLZ2(x,f,nvar,nobj);

		double r = 0.5;
		if (nobj == 2) r = 0.15;
		else if (nobj == 3) r = 0.4;

		double max_i =  1<<30;
		double tmp1 = 0;
		double tmp2 = 0;
		double sum = 0;
		for (int i = 0; i < nobj; i++)
		{
			sum += pow(f[i],2.0);
		}

		for (int i = 0; i < nobj; i++)
		{
			tmp1 = pow(f[i] - 1, 2.0) + sum - pow(f[i], 2.0) - pow(r, 2.0);
			if (tmp1 < max_i)
				max_i = tmp1;

			tmp2 += pow(f[i] - 1.0 / sqrt((double)nobj), 2.0);
		}

		tmp2 -= pow(r, 2.0);

		if (tmp2 < max_i)
			c[0] = -tmp2;
		else
			c[0] = -max_i;
	}

	void C2_Convex_DTLZ2(double *x, double *f, double *c, const unsigned int nvar, const unsigned int nobj)
	{
		DTLZ::Convex_DTLZ2(x, f, nvar, nobj);

		double r[] = { 0.1, 0.225, 0.225, 0.225, 0.225, 0.225, 0.26, 0.26, 0.26, 0.26, 0.26, 0.26, 0.26, 0.27 };

		double sum = 0;
		for (int i = 0; i < nobj; i++)
		{
			sum += f[i];
		}
		double lamda = 1.0 / nobj * sum;
		double tmp = 0;
		for (int i = 0; i < nobj; i++)
		{
			tmp += pow(f[i]-lamda,2.0);
		}

		c[0] = tmp - pow(r[nobj - 2], 2.0);
	}

	void C3_DTLZ1(double *x, double *f, double *c, const unsigned int nvar, const unsigned int nobj)
	{
		DTLZ::DTLZ1(x, f, nvar, nobj);

		double sum = 0;
		for (int i = 0; i < nobj; i++)
		{
			sum += 2*f[i];
		}

		for (int i = 0; i < nobj; i++)
		{
			c[i] = sum - f[i] - 1;
		}
	}
	void C3_DTLZ4(double *x, double *f, double *c, const unsigned int nvar, const unsigned int nobj)
	{
		DTLZ::DTLZ4(x, f, nvar, nobj);

		double sum = 0;
		for (int i = 0; i < nobj; i++)
		{
			sum += pow(f[i], 2.0);
		}

		for (int i = 0; i < nobj; i++)
		{
			c[i] = sum - 0.75 * pow(f[i], 2.0) - 1;
		}
	}
};