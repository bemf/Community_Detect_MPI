
#ifndef __DHEAInd_CLASS_H_
#define __DHEAInd_CLASS_H_

#include "../global.h"


class CHEAInd{
public:
    //vector <double> x_var;
	vector <double> x_var;
	vector <double> y_obj;

	int sectorialindex;
	double k_value;

    bool   record_update;
    void   rnd_init();
    void   obj_eval();
	void   obj_eval_graph(MGraph *G);
    void   show_objective();

	//图
	void   rnd_init(MGraph *G);

    CHEAInd();
	~CHEAInd();

	void cal_k_value(const vector<double>& idealPoint,const int& hyperplane_intercept);
	vector<double> cal_normalizedf(const vector<double>& idealPoint,const int& hyperplane_intercept);
	void cal_V_obj(const vector<double>& idealPoint,const int& hyperplane_intercept, vector<int> &V_obj);

	//计算个体的对应子问题编号
	inline void obj_index(const vector <double> & idealPoint,const int& hyperplane_intercept)
	{
		vector<int> V_obj(nobj,0);
		cal_V_obj(idealPoint, hyperplane_intercept, V_obj);
		sectorialindex = cal_index_from_lambda(V_obj, hyperplane_intercept);//get_index_from_V_obj(V_obj, hyperplane_intercept);
        //if(sectorialindex>=100)
        //{
            //show_objective();
            //std::cout<<"IdealPoint:["<<idealPoint[0]<<","<<idealPoint[1]<<"] ";
            //std::cout<<"hyperplane_Intercept:"<<hyperplane_intercept<<" ";
            //std::cout<<"V_obj:["<<V_obj[0]<<","<<V_obj[1]<<"] ";
            //std::cout<<"sectorialindex:"<<sectorialindex<<endl;
        //}
		V_obj.clear();
	}
	//通过观察向量计算对应子问题的编号
	int get_index_from_V_obj(const vector<int>& V_obj, const int& hyper_Intercept);


	/** Compares two points in multiple objective space
	*	Returns _Dominated if this is dominated
	*	Returns _Dominating if this is dominating
	*	It may also return _Nondominated or _Equal */

	TCompare Compare(CHEAInd& ind2);

	void   operator=(const CHEAInd &ind2);
};

//子问题定义
class CHEA_SOP{
public:
	vector <int> V_obj_o;//以Intercept为截距的超平面上的中心观察向量
	int sectorialindex;
	CHEAInd indiv;

	vector<int> neighborhood;

	//获取子问题的邻居子问题编号
	/*
	inline void get_vicinity(int vicinity_range, int hyperplane_Intercept)
	{
		vector<int> vicinity(nobj, 0);
		cal_vicinity(V_obj_o, vicinity_range, hyperplane_Intercept, 0, 0, false, vicinity, neighborhood);
		vicinity.clear();

	}
	*/
	inline void get_vicinity(int hyperplane_Intercept)
	{
		//vector<int> vicinity(nobj, 0);
		//cal_vicinity(V_obj_o, vicinity_range, hyperplane_Intercept, 0, 0, false, vicinity, neighborhood);
		//vicinity.clear();
		vector<int> vicinity= V_obj_o;
		for(int i=0;i<nobj;i++)
		{
			if(vicinity[i]==0)
				continue;
			vicinity[i] -= 1;
			for(int j=0;j<nobj;j++)
			{
				if(j == i)
					continue;
				if(vicinity[j]== hyperplane_Intercept)
					continue;
				vicinity[j] += 1;

				int index = indiv.get_index_from_V_obj(vicinity, hyperplane_Intercept);
				neighborhood.push_back(index);
				vicinity[j] -= 1;
			}
			vicinity[i] += 1;
		}
	}

	//计算邻居子问题
	void cal_vicinity(vector<int> &V_obj, int vicinity_range,int hyperplane_Intercept, int cal_index_now, int left_range, bool is_change_before, vector<int> vicinity, vector<int> &vicinity_set);


	CHEA_SOP();
	~CHEA_SOP();

};

#endif
