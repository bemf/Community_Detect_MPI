#include "common.h"


void minfastsort(vector<double> &x, vector<int> &idx, int n, int m)
{
	minfastsort((double*)&(*(x.begin())), (int*)&(*(idx.begin())),  n,  m);
}

void minfastsort(double* x, int* idx, int n, int m)
{
	int selectIndex;
	for ( int i = 0; i < m; i++)
	{
		selectIndex = i;
		for ( int j = i + 1; j < n; j++)
		{
			if (x[j] < x[selectIndex])
			{
				selectIndex = j;
			}
		}
		if (i != selectIndex)
		{
			double temp = x[selectIndex];
			x[selectIndex] = x[i];
			x[i] = temp;
			temp = idx[selectIndex];
			idx[selectIndex] = idx[i];
			idx[i] = temp;
		}
	}
}


double dist_vector(vector <double> &vec1, vector <double> &vec2)
{
	int dim = vec1.size();
    double sum = 0;
	for(int n=0; n<dim; n++)
	    sum+=(vec1[n] - vec2[n])*(vec1[n] - vec2[n]);
	return sqrt(sum);
}

double observepoint_dist_vector(vector <double> &vec1, vector <double> &vec2)
{
	int dim = vec1.size();
	double sum = 0;
	for(int n=0; n<dim; n++)
		sum+=(vec1[n] - (vec2[n] ) )*(vec1[n] - (vec2[n] ));
	return sqrt(sum);
}

//这个有错
void quicksort_increasing(vector <vector<double> > &paretofront, int nPoints)
{
	double temp;
	for (int i = 0; i<nPoints; i++)
	{
		for (int j = i + 1; j<nPoints; j++)

		if (paretofront[i][0]>paretofront[j][0])
		{
			temp = paretofront[i][0];
			paretofront[i][0] = paretofront[j][0];
			paretofront[j][0] = temp;

			temp = paretofront[i][1];
			paretofront[i][1] = paretofront[j][1];
			paretofront[j][1] = temp;
		}
	}
}


TCompare ParetoCompare(vector <double> & y_obj1, vector <double>& y_obj2) {
	bool bBetter = false;
	bool bWorse = false;

	int i = 0;
	do {
		if (y_obj1[i] < y_obj2[i])
			bBetter = true;
		if (y_obj2[i] < y_obj1[i])
			bWorse = true;
		i++;
	} while (!(bWorse && bBetter) && (i < nobj));
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

void save_pop_front(vector <vector<double> >  &pop_front, char saveFilename[1024])
{
	std::fstream fout;
	fout.open(saveFilename, std::ios::out);
	for (int n = 0; n<pop_front.size(); n++)
	{
		for (int k = 0; k<nobj; k++)
		{
			fout << pop_front[n][k] << "  ";
		}
		fout << "\n";
	}
	fout.close();
}


void random_permutation(int *perm, int size)
{
	int *index = new int[size];
	bool *flag = new bool[size];
	for (int n = 0; n<size; n++)  {
		index[n] = n;
		flag[n] = true;
	}

	int num = 0;
	while (num<size){
		int start = int(size*rnd_uni(&rnd_uni_init));
		while (1){
			if (flag[start]){
				perm[num] = index[start];
				flag[start] = false;
				num++;
				break;
			}
			if (start == (size - 1))
				start = 0;
			else
				start++;
		}
	}

	delete[] index;
	delete[] flag;

}

//计算组合数
int choose(int n, int m)
{
	/*
	int numerator = 1;
	int denominator = 1;
	for (int i = n; i > n - m; i--)
		numerator *= i;
	for (int i = m; i > 1; i--)
		denominator *= i;
	return numerator / denominator;
	*/
	if (m > n)
    {
        return 0;
    }
    return (int)floor(0.5 + exp(lnchoose(n, m)));
}

double lnchoose(int n, int m)
{
    if (m > n)
    {
        return 0;
    }
    if (m < n/2.0)
    {
        m = n-m;
    }
    double s1 = 0;
    for (int i=m+1; i<=n; i++)
    {
        s1 += log((double)i);
    }
    double s2 = 0;
    int ub = n-m;
    for (int i=2; i<=ub; i++)
    {
        s2 += log((double)i);
    }
    return s1-s2;
}


int cal_index_from_lambda(const vector<int> &V_obj,int n_H)
{
	vector<int> h(nobj - 1, 0);
	h[nobj - 2] = n_H - V_obj[nobj - 1];
	for (int i = nobj - 3; i >0; i--)
	{
		h[i] = h[i + 1] - V_obj[i + 1];
	}
	h[0] = V_obj[0];

	int result_index = 0;
	for (int i = 0; i < nobj - 1; i++)
	{
		result_index += choose(h[i] + i, i + 1);
	}
	h.clear();

	return result_index;
}

vector<int> gen_vicinities(vector<int> V_obj,int n_H)
{
	/*
	vector<int> vicinities;
	for (int i = 0; i < nobj; i++)
	{
	V_obj[i] = V_obj[i] + 1;
	if (V_obj[i] <= n_H)
	{
	for (int j = i + 1; j < nobj; j++)
	{
	V_obj[j] = V_obj[j] - 1;
	if (V_obj[j] >= 0)
	{
	int index = cal_index_from_lambda(V_obj);
	vicinities.push_back(index);
	}
	V_obj[j] = V_obj[j] + 1;
	}
	}
	V_obj[i] = V_obj[i] - 2;
	if (V_obj[i] >= 0)
	{
	for (unsigned int j = i + 1; j < nobj; j++)
	{
	V_obj[j] = V_obj[j] + 1;
	if (V_obj[j] <= n_H)
	{
	int index = cal_index_from_lambda(V_obj);
	vicinities.push_back(index);
	}
	V_obj[j] = V_obj[j] - 1;
	}
	}
	V_obj[i] = V_obj[i] + 1;
	}
	*/

	vector<int> vicinities;
	for (int i = 0; i < nobj; i++)
	{
		V_obj[i] = V_obj[i] + 1;
		if (V_obj[i] <= n_H)
		{
			for (int j = 0; j < nobj; j++)
			{
				if (j == i) continue;

				V_obj[j] = V_obj[j] - 1;
				if (V_obj[j] >= 0)
				{
					int index = cal_index_from_lambda(V_obj,n_H);
					vicinities.push_back(index);
				}
				V_obj[j] = V_obj[j] + 1;
			}
		}
		V_obj[i] = V_obj[i] - 1;
	}

	return vicinities;
}


void population2front(vector<vector<double> >&pop_y, vector<vector<double> >&pf_y)
{
	vector<int> nDominated;
	for (int n = 0; n < pop_y.size(); n++)
		nDominated.push_back(0);

	for (int k = 0; k<pop_y.size(); k++)
		for (int j = k + 1; j<pop_y.size(); j++)
		{
			TCompare tresult = ParetoCompare(pop_y[k], pop_y[j]);
			if (tresult == _Pareto_Dominated)
				nDominated[k]++;
			else if (tresult == _Pareto_Dominating)
				nDominated[j]++;
		}

	for (int n = 0; n<pop_y.size(); n++)
		if (nDominated[n] == 0)
			pf_y.push_back(pop_y[n]);
	nDominated.clear();
}

void population2front(vector<vector<double> >&pop_y, vector<vector<double> >&pf_y, vector<int> &pop_index)
{
	vector<int> nDominated;
	for (int n = 0; n<pop_y.size(); n++)
		nDominated.push_back(0);


	for (int k = 0; k<pop_y.size(); k++)
		for (int j = k + 1; j<pop_y.size(); j++)
		{
			TCompare tresult = ParetoCompare(pop_y[k], pop_y[j]);
			if (tresult == _Pareto_Dominated)
				nDominated[k]++;
			else if (tresult == _Pareto_Dominating)
				nDominated[j]++;
		}

	for (int n = 0; n<pop_y.size(); n++)
		if (nDominated[n] == 0) {
			pf_y.push_back(pop_y[n]);
			pop_index.push_back(n);
		}
	nDominated.clear();
}

void save_population_front(vector <vector<double> >  &pf_y, char saveFilename[1024])
{
	std::fstream fout;
	fout.open(saveFilename, std::ios::out);
	for (int n = 0; n<pf_y.size(); n++)
	{
		for (int k = 0; k<nobj; k++)
		{
			fout << pf_y[n][k] << "  ";
		}
		fout << "\n";
	}
	fout.close();
}

void save_population_front_split(vector <vector<double> >  &pf_y, vector<int>&pf_index, vector<int>&index_split, char saveFilename[1024])
{
	std::fstream fout;
	fout.open(saveFilename, std::ios::out);
	for (int n = 0; n<pf_y.size(); n++)
	{
		fout << index_split[pf_index[n]] << "  ";
		for (int k = 0; k<nobj; k++)
		{
			fout << pf_y[n][k] << "  ";
		}
		fout << "\n";
	}
	fout.close();
}


void save_pos_split(vector <vector<double> >  &pos, vector<int>&pf_index, vector<int>&index_split, char saveFilename[1024])
{
	std::fstream fout;
	fout.open(saveFilename, std::ios::out);
	for (int n = 0; n<pf_index.size(); n++)
	{
		fout << index_split[pf_index[n]] << "  ";
		for (int k = 0; k<nvar; k++)
			fout << pos[n][k] << "  ";
		fout << "\n";
	}
	fout.close();
}



