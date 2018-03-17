#include "../common.h"
//#include "cec09.h"
#include "cheaind.h"
#include "../TestFunction/testInstance.h"

CHEAInd::CHEAInd()
{

    for(int i=0; i<nvar; i++)
        x_var.push_back(0.0);
    for(int j=0; j<nobj; j++)
        y_obj.push_back(0.0);


    record_update=false;
}

CHEAInd::~CHEAInd()
{
    x_var.clear();
    y_obj.clear();
}

//图初始化
void CHEAInd::rnd_init(MGraph *G)
{
	//for (int i = 0; i < G->numVertexes; i++){
		//int temp = 0;

		//int rank = int(rnd_uni(&rnd_uni_init) * numA[i]) + 1;
		//for (int j = 0; j < G->numVertexes; j++){
			//if (G->arc[i][j] == 1){
				//temp++;
				//if (rank == temp){
					//x_var[i]=j;
				//}
			//}
		//}
	//}
	for (int i = 0; i < G->numVertexes; i++)
    {
        if(G->degree[i]>0)
        {
    		int j= int(rnd_uni(&rnd_uni_init) * G->degree[i]) ;
            //cout <<"degree "<<G->degree[i]<<endl;;
            //cout <<"i "<<i<<" j "<<j<<endl;
            //cout<<" nodes "<<G->Adjacency[i][j]<<endl;
    		x_var[i]=G->Adjacency[i][j];

        }
        else
    		x_var[i]=i;
        //cout<<"i "<<i<<" x_var[i] "<<x_var[i]<<endl;
    }

	obj_eval_graph(G);
}

void CHEAInd::obj_eval_graph(MGraph *G)
{
	objectives_graph2(x_var, y_obj,G);
}

//个体初始化
void CHEAInd::rnd_init()
{
    for(int n=0;n<nvar;n++)
	{
        x_var[n] = lowBound[n] + rnd_uni(&rnd_uni_init)*(uppBound[n] - lowBound[n]);
	}

	obj_eval();
}

//获取目标值
void CHEAInd::obj_eval()
{
	objectives(x_var, y_obj);
}

void CHEAInd::show_objective()
{
	for(int n=0;n<nobj;n++)
		std::cout<<y_obj[n]<<" ";
		std::cout<<sectorialindex<<" ";
	std::cout<<"\n";

}

void CHEAInd::operator=(const CHEAInd &ind2)
{
    x_var = ind2.x_var;
    y_obj = ind2.y_obj;


    sectorialindex=ind2.sectorialindex;
    k_value=ind2.k_value;

    record_update=ind2.record_update;

}

//支配关系比较
TCompare CHEAInd::Compare(CHEAInd& ind2) {
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

void CHEAInd::cal_k_value(const vector<double>& idealPoint,const int& hyperplane_intercept)
{
	k_value = 0.0;
	//标准化
	for (int j = 0; j<nobj; j++)
	{
		k_value += (y_obj[j] - idealPoint[j]);
	}

	if (k_value == 0 )
		k_value = hyperplane_intercept;
}

vector<double> CHEAInd::cal_normalizedf(const vector<double>& idealPoint,const int& hyperplane_intercept)
{
	double obj_relative_sum = 0.0;
	vector<double> normalizedf(nobj,0);
	//normalizedf.reserve(nobj);
	//标准化
	for (int j = 0; j<nobj; j++)
	{
		normalizedf[j] = y_obj[j] - idealPoint[j];
		obj_relative_sum += normalizedf[j];
	}

	if (obj_relative_sum == 0 )
	{
		obj_relative_sum = hyperplane_intercept;
		normalizedf[0] = hyperplane_intercept;
		for (int j =1; j < nobj; j++)
		{
			normalizedf[j] = 0;
		}
	}
	else
	{//映射到超平面上的点坐标
		for (int j = 0; j < nobj; j++)
		{
			normalizedf[j] = normalizedf[j] / obj_relative_sum * hyperplane_intercept;
		}
	}
	k_value = obj_relative_sum;
	return normalizedf;
}

//计算目标点对应子问题的观察向量
void CHEAInd::cal_V_obj(const vector<double>& idealPoint,const int& hyperplane_intercept, vector<int> &V_obj_o)
{
	//vector<double> normalizedf = cal_normalizedf(idealPoint,hyperplane_intercept);
	double obj_relative_sum = 0.0;
	vector<double> normalizedf(nobj,0);
	//normalizedf.reserve(nobj);
	//标准化
	for (int j = 0; j<nobj; j++)
	{
		normalizedf[j] = y_obj[j] - idealPoint[j];
		obj_relative_sum += normalizedf[j];
	}

	if (obj_relative_sum == 0 )
	{
		obj_relative_sum = hyperplane_intercept;
		normalizedf[0] = hyperplane_intercept;
		for (int j =1; j < nobj; j++)
		{
			normalizedf[j] = 0;
		}
	}
	else
	{//映射到超平面上的点坐标
		for (int j = 0; j < nobj; j++)
		{
			normalizedf[j] = normalizedf[j] / obj_relative_sum * hyperplane_intercept;
		}
	}
	k_value = obj_relative_sum;

	vector<bool> is_complete_bit(nobj, false);
	int size = nobj;
	int min_index = 0;
	double min_value = 1e7;
	vector<double> lowdif(nobj,0);
	vector<double> uppdif(nobj,0);
	double fix;
	double selectdata;
	for (int i = 0; i < nobj; i++)
	{
		lowdif[i] =	(double)floor(normalizedf[i])- normalizedf[i];
		uppdif[i] = 1.0 + lowdif[i];
		selectdata = (-lowdif[i]) < uppdif[i] ? - lowdif[i] : uppdif[i];
		if (min_value > selectdata)
		{
			min_index = i;
			min_value = selectdata;
		}
	}

	while (1)
	{
		selectdata = ((-lowdif[min_index]) < uppdif[min_index] ? lowdif[min_index] : uppdif[min_index]);
		V_obj_o[min_index] = selectdata<0 ? floor(normalizedf[min_index]) : ceil(normalizedf[min_index]);

		is_complete_bit[min_index] = true;
		size--;

		if (size == 1)
		{
			int sum = 0;
			int left_index;
			for (int i = 0; i < nobj; i++)
			{
				if (is_complete_bit[i])
					sum += V_obj_o[i];
				else
					left_index = i;
			}
			V_obj_o[left_index] = hyperplane_intercept - sum;
			break;
		}
		fix =  selectdata / size;

		min_value = 1e7;
		for (int i = 0; i < nobj; i++)
		{
			if (is_complete_bit[i])continue;

			lowdif[i] += fix;
			uppdif[i] += fix;
			selectdata = (-lowdif[i]) < uppdif[i] ? -lowdif[i] : uppdif[i];
			if (min_value > selectdata)
			{
				min_index = i;
				min_value = selectdata;
			}
		}
	}

	is_complete_bit.clear();
	lowdif.clear();
	uppdif.clear();

	normalizedf.clear();
}


//由观察向量坐标获取对应子问题的编号
int CHEAInd::get_index_from_V_obj(const vector<int> &V_obj,const int& hyper_Intercept)
{
	/*
	if (2 == nobj) return V_obj[0];//二维目标
	else if (3 == nobj) return (hyper_Intercept - V_obj[2] + 1)*(hyper_Intercept - V_obj[2]) / 2 + V_obj[0];//三维目标

	vector<int> h(nobj-1,0);
	h[0] = hyper_Intercept - V_obj[nobj-1];
	for (int i = 1; i < nobj-1 ; i++)
	{
		h[i] = h[i - 1] - V_obj[nobj-1-i];
	}

	int result_index = 0;
	for (int i = 0; i < nobj - 1; i++)
	{
		result_index += choose(h[i] + nobj -2 -i,nobj-1-i);
	}
	h.clear();
	return result_index;
	*/
	if (2 == nobj) return V_obj[0];//二维目标
	else if (3 == nobj) return (hyper_Intercept - V_obj[2] + 1)*(hyper_Intercept - V_obj[2]) / 2 + V_obj[0];//三维目标

	vector<int> h(nobj-1,0);
	h[nobj-2] = hyper_Intercept - V_obj[nobj-1];
	for (int i = nobj-3; i >0; i--)
	{
		h[i] = h[i+1] - V_obj[i+1];
	}
	h[0] = V_obj[0];

	int result_index = 0;
	for (int i = 0; i < nobj - 1; i++)
	{
		result_index += choose(h[i]+i ,i+1);
	}
	h.clear();
	return result_index;
	/*
	if (2 == nobj) return V_obj[0];//二维目标
	else if (3 == nobj) return (hyper_Intercept - V_obj[2] + 1)*(hyper_Intercept - V_obj[2]) / 2 + V_obj[0];//三维目标
	else return cal_index(V_obj,nobj-1,hyper_Intercept,true) -1;//更高维度目标
	*/
}

CHEA_SOP::CHEA_SOP()
{
	V_obj_o = vector<int>(nobj,0);
}
CHEA_SOP::~CHEA_SOP()
{
	V_obj_o.clear();
}

//计算子问题的邻居子问题编号
void CHEA_SOP::cal_vicinity(vector<int> &V_obj, int vicinity_range, int hyperplane_Intercept,int cal_index_now, int left_range, bool is_change_before, vector<int> vicinity, vector<int> &vicinity_index_set)
{
	int index_value = V_obj[cal_index_now];
	for (int k = index_value - vicinity_range; k <= index_value + vicinity_range; k++)
	{
		if (k <0 || k > hyperplane_Intercept) continue;
		vicinity[cal_index_now] = k;
		if (cal_index_now == nobj - 1)
		{
			if (is_change_before && (left_range + (k - index_value) == 0))
			{
				int index = indiv.get_index_from_V_obj(vicinity, hyperplane_Intercept);
				vicinity_index_set.push_back(index);
			}
			continue;
		}
		bool true_or_false = true;
		if (!is_change_before && (k == index_value)) true_or_false = false;
		cal_vicinity(V_obj, vicinity_range, hyperplane_Intercept, cal_index_now + 1, left_range + (k - index_value), true_or_false, vicinity, vicinity_index_set);
	}
}
