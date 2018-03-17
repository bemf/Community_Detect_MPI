#include "ccaeaalg.h"


bool CAEA::Pareto_HyperVolume_compare_sectorialgrid(const CHEAInd& ind, CHEAInd& replacedInd)
{
	bool bReplaced = false;

	double contributing1, contributing2;
	CHEA_SOP &subproblem = population[ind.sectorialindex];
	subproblem.indiv.obj_index(IdealPoint, hyperplane_Intercept);
	if (subproblem.sectorialindex == subproblem.indiv.sectorialindex)
	{
		//在此处设置动态参考点或者静态参考点

		contributing1 = GetHyperVolume(ind, ReferencePoint);//
		contributing2 = GetHyperVolume(subproblem.indiv, ReferencePoint);//

		if (contributing1 < contributing2)
		{
			replacedInd = population[ind.sectorialindex].indiv;
			population[ind.sectorialindex].indiv = ind;//replace
			bReplaced = true;
		}
	}
	else {
		//sectorpop[ind.sectorialindex].indiv = ind;//replace
		//bReplaced = true;
		replacedInd = population[ind.sectorialindex].indiv;
		population[ind.sectorialindex].indiv = ind;
		//Pareto_HyperVolume_compare_sectorialgrid(replacedInd);// , replacedInd);
		bReplaced = true;
	}

	return bReplaced;
}

double CAEA::GetHyperVolume(const CHEAInd&  ind, vector <double> &ref_ponit_cal)
{
	double FastigiateVolume;
	int ind_index = ind.sectorialindex;

	double normalizedf[2];
	for (int j = 0; j<nobj; j++)
	{
		normalizedf[j] = ind.y_obj[j] - IdealPoint[j];
	}
	if (ind_index == popsize -1) {
		FastigiateVolume = 0.5*(ind_index + 0.5) / (popsize - ind_index - 1.5)*normalizedf[1] * normalizedf[1] + (ReferencePoint[1] - normalizedf[1])*normalizedf[0];
	}
	else if (ind_index == 0) {
		FastigiateVolume = 0.5*(popsize - ind_index - 0.5) / (ind_index - 0.5)*normalizedf[0] * normalizedf[0] + (ReferencePoint[0] - normalizedf[0])*normalizedf[1];
	}
	else {
		FastigiateVolume = 0.5*(popsize - ind_index - 0.5) / (ind_index - 0.5)*normalizedf[0] * normalizedf[0] + 0.5*(ind_index + 0.5) / (popsize - ind_index - 1.5)*normalizedf[1] * normalizedf[1] - normalizedf[0] * normalizedf[1];
	}

	return FastigiateVolume;
}

