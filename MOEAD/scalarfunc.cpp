#include "scalarfunc.h"

double norm_vector(vector <double> &x)
{
	double sum = 0;
	for (int i = 0; i<x.size(); i++)
		sum = sum + x[i] * x[i];
	return sqrt(sum);
}

double innerproduct(vector <double>&vec1, vector <double>&vec2)
{
	double sum = 0;
	for (int i = 0; i<vec1.size(); i++)
		sum += vec1[i] * vec2[i];
	return sum;
}


// scalarizing functions for decomposition methods
double scalar_func(vector <double> &y_obj, vector <double> &namda, vector<double> &idealpoint, CIndividual* nbi_node)
{

	double fvalue = 0;

	// Tchebycheff approach
	if (tDecompType == _Decomposition_TCH1)
	{
		double max_fun = -1.0e+30;
		for (int n = 0; n < nobj; n++)
		{
			double diff = fabs(y_obj[n] - idealpoint[n]);

			double feval;
			if (namda[n] == 0)
				feval = 0.00001*diff;
			else
				feval = diff*namda[n];
			if (feval>max_fun) max_fun = feval;

		}

		fvalue = max_fun;
	}

	// normalized Tchebycheff approach
	if (tDecompType == _Decomposition_TCH2)
	{
		vector <double> scale;
		for (int i = 0; i<nobj; i++)
		{
			double min = 1.0e+30, max = -1.0e+30;
			for (int j = 0; j<nobj; j++)
			{
				double tp = nbi_node[j].y_obj[i];
				if (tp>max) max = tp;
				if (tp<min) min = tp;
			}
			scale.push_back(max - min);
			if (max - min == 0) return 1.0e+30;
		}

		double max_fun = -1.0e+30;
		for (int n = 0; n<nobj; n++)
		{
			double diff = (y_obj[n] - idealpoint[n]) / scale[n];
			double feval;
			if (namda[n] == 0)
				feval = 0.0001*diff;
			else
				feval = diff*namda[n];
			if (feval>max_fun) max_fun = feval;

		}
		fvalue = max_fun;
	}


	//* Boundary intersection approach
	if (tDecompType == _Decomposition_PBI)
	{

		// normalize the weight vector (line segment)
		double nd = norm_vector(namda);
		///for(int i=0; i<nobj; i++)
		///	namda[i] = namda[i]/nd;

		vector <double> realA(nobj);
		vector <double> realB(nobj);

		// difference beween current point and reference point
		for (int n = 0; n<nobj; n++)
			realA[n] = (y_obj[n] - idealpoint[n]);

		// distance along the line segment
		double d1 = fabs(innerproduct(realA, namda)) / nd;

		// distance to the line segment
		for (int n = 0; n<nobj; n++)
			realB[n] = (y_obj[n] - (idealpoint[n] + d1*namda[n]));
		double d2 = norm_vector(realB);

		fvalue = d1 + 5 * d2;

	}

	return fvalue;
}