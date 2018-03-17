#pragma once
#include "../CHEA/cheaalg.h"

class CAEA : public CHEA
{
public:
	CAEA(int pop_size, int hyper_intercept) :CHEA(pop_size, hyper_intercept) {}

	bool Pareto_HyperVolume_compare_sectorialgrid(const CHEAInd& ind, CHEAInd& replacedInd);

	double GetHyperVolume(const CHEAInd&  ind, vector <double> &ref_ponit_cal);

};
