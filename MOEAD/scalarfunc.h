#ifndef __SCALARFUNC_H_
#define __SCALARFUNC_H_

#include "../global.h"
#include "moeadind.h"

double norm_vector(vector <double> &x);

double innerproduct(vector <double>&vec1, vector <double>&vec2);

// scalarizing functions for decomposition methods
double scalar_func(vector <double> &y_obj, vector <double> &namda, vector<double> &idealpoint, CIndividual* nbi_node);


#endif
