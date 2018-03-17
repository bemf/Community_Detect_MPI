
#ifndef __CAEA_CLASS_H_
#define __CAEA_CLASS_H_

#include "../global.h"

class CCAEAInd{
public:
    vector <double> x_var;
	vector <double> y_obj;

	vector <int> x_graph;
	vector <double> y_graph;
    //int    rank, count;
	//vector <int>        sectorialcoord;
	//vector <double>     
	double sectorialangle;
	int sectorialindex;

    void   rnd_init();
	void   rnd_init(MGraph *G);
    void   obj_eval();
	void   obj_eval(MGraph *G);
    void   show_objective();
	CCAEAInd();
    CCAEAInd(MGraph *G);
	~CCAEAInd();



	//void  yobj2angle(vector <double> & observerpoint);

	void obj2angle_index(vector<double> & IdealPoint, int sectornum);

	//bool angle_in_sector(vector <double> & observerpoint, int sectornum, double anglesingle ,int sectorcoord);


		/** Compares two points in multiple objective space
	*
	*	Returns _Dominated if this is dominated
	*	Returns _Dominating if this is dominating
	*	It may also return _Nondominated or _Equal */
	TCompare Compare(CCAEAInd& ind2);

	//bool   operator<(const CCAEAInd& ind2);
	//bool   operator==(const CCAEAInd& ind2);
	//void   operator=(const CCAEAInd& ind2);

};

#endif