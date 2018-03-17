//#ifndef __GLOBAL_H_
//#include "global.h"
//#endif
//#include "cec09.h"

#include "moeadind.h"
#include "../TestFunction/testInstance.h"

CIndividual::CIndividual()
{
	x_var = vector<double>(nvar, 0);
    y_obj = vector<double>(nobj, 0);
}

CIndividual::~CIndividual()
{
	x_var.clear();
	y_obj.clear();
}

void CIndividual::rnd_init()
{
    for(int n=0;n<nvar;n++)
        x_var[n] = lowBound[n] + rnd_uni(&rnd_uni_init)*(uppBound[n] - lowBound[n]);   

	obj_eval();//后来加的
}

void CIndividual::obj_eval()
{
	//obj_eval_functions(x_var, y_obj);
	objectives(x_var, y_obj);	
}


void CIndividual::show_objective()
{
    for(int n=0; n<nobj; n++)
		printf("%f ",y_obj[n]);
	printf("\n");
}

void CIndividual::show_variable()
{
    for(int n=0; n<nvar; n++)
		printf("%f ",x_var[n]);
	printf("\n");
}

void CIndividual::operator=(const CIndividual &ind2)
{
    x_var = ind2.x_var;
	y_obj = ind2.y_obj;
}

bool CIndividual::operator<(const CIndividual &ind2)
{
	bool dominated = true;
    for(int n=0; n<nobj; n++)
	{
		if(ind2.y_obj[n]<y_obj[n]) return false;
	}
	if(ind2.y_obj==y_obj) return false;
	return dominated;
}


bool CIndividual::operator<<(const CIndividual &ind2)
{
	bool dominated = true;
    for(int n=0; n<nobj; n++)
	{
		if(ind2.y_obj[n]<y_obj[n]  - 0.0001) return false;
	}
	if(ind2.y_obj==y_obj) return false;
	return dominated;
}

bool CIndividual::operator==(const CIndividual &ind2)
{
	if(ind2.y_obj==y_obj) return true;
	else return false;
}

TCompare CIndividual::Compare(CIndividual& ind2) {
	bool bBetter = false;
	bool bWorse = false;

	int i = 0;
	do { 
		if(y_obj[i] < ind2.y_obj[i])
			bBetter = true;
		if(ind2.y_obj[i] < y_obj[i])
			bWorse = true;
		i++;
	}while (!(bWorse && bBetter) && (i < nobj));
	if (bWorse) {
		if (bBetter)
			return _Pareto_Nondominated;
		else
			return _Pareto_Dominated;
	}
	else {
		if (bBetter)
			return _Pareto_Dominating;
		else
			return _Pareto_Equal;
	}
}


// defining subproblem 


TSOP::TSOP()
{
}

TSOP::~TSOP()
{
}


void TSOP::operator=(const TSOP&sub2)
{
	indiv  = sub2.indiv;
	table  = sub2.table;
	namda  = sub2.namda;
	array  = sub2.array;
}



//#endif