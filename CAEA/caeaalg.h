#ifndef __CAEA_FUNC_H_
#define __CAEA_FUNC_H_


#include "../global.h"
#include "../common.h"
#include "../recomb.h"

#include "../wfg.h"

#include "caeaind.h"

class CCAEA
{
public:
	CCAEA(int pop_size);
	virtual ~CCAEA();

	void initTopologies();
	void execute(int mg, int irun, vector<double>& igd, vector<double>& hv, double &runtime, PTYPE ptype);
	void execute(int mg, int irun, vector<double>& igd, vector<double>& hv, double &runtime, PTYPE ptype,MGraph *G);
	void init_population();
	void init_population(MGraph *G);

	void save_front(char savefilename[1024]);          // save the pareto front into files
	void save_pos(char saveFilename[1024]);

	//vector <int>  sectorialgrid;

	vector <CCAEAInd>  population;

	vector <CCAEAInd>  ps;

	CCAEAInd onechild;

	//int sectornum;
	int popsize;
	//double anglesingle;

	//int p_nparts;

	//vector<vector <double>> AnchorPoint;
	vector <double> TrueNadirPoint;

    //vector <double> UtopiaPoint;
	//vector <double> PseudoNadirPoint;

	vector <double> IdealPoint;

	vector <double> ReferencePoint;

	//double distance;

	//int    nfes;

    //void operator=(const CCAEA &alg);

	bool Pareto_HyperVolume_compare_sectorialgrid(CCAEAInd& ind);

	//void initial_observation_and_reference_point(CCAEAInd& ind);
	bool update_extreme_point(CCAEAInd& ind);
	bool update_ideal_point(const double* y_obj);

	//void uniform_selection(CCAEAInd*& ind_selected);

	void reset_angle();

	//void sbx_mutation_evolution(CCAEAInd& ind_selectedparents1, CCAEAInd& ind_selectedparents2);

	double GetFastigiateHyperVolume(CCAEAInd&  ind, int ind_index,vector<double> &ReferentPoint);

	int  tour_selection_hv2(vector <CCAEAInd>  &mypopulation);
	double  tour_selection_hv_difference(int p, vector <CCAEAInd>  &mypopulation);


	void operator=(const CCAEA &emo);

	int *rcounts_y;
	int* displs_y;
	int *rcounts_a;
	int* displs_a;

	void gather_pop_y();
	void gather_populations();
	void calc_qualities(double& igd, double &hv);
	double calc_distance(vector<vector<double> >& ParetoFront);
	double calc_hypervolume(vector<vector<double> >& ParetoFront);

	clock_t start, finish;
	clock_t unevolvetime;


	bool update_idealpoint_from_outside();
	void sendIdealpoint(double *buf, bool &isSend, int target, MPI_Request* request);
	void getAndUpdateIndivs(double* b_recv_indiv, int datasize, int source, MPI_Request *request);
	void update_indiv_from_outside();

	//Parallel Parameter
	vector<int>neighbor_partitions;
	int partition_niche;

	//int *pops_sub;
	int *pops_sub_self;
	int *pops_sub_start;
	double *b_recv, *b_send;
	int s_recv, s_send;

	//更新理想点
	bool isIdealpointUpdate;
	MPI_Request *req_idealpoint_recv;
	MPI_Request *req_idealpoint_send;
	bool *isBegin_Ideal;
	double **recv_idealpoint, **send_idealpoint;
	MPI_Status status_idealpoint;

	//更新个体
	bool *isBegin_Indiv;
	bool *recordUpdateofIndiv;
	MPI_Request *req_indiv_recv;
	MPI_Request *req_indiv_send;
	double **recv_indiv, **send_indiv;
	double *indiv_buf;
	int s_indiv;
	int indivBufSize;
	MPI_Status status_indiv;

	//完成状态
	MPI_Request *req_done_recv;
	MPI_Request *req_done_send;
	int *isDoneRecv;
	int *isDoneSend;
	int *b_isdone;
	MPI_Status status_done;

};

void population2front(vector <CCAEAInd>  &mypopulation, vector <CCAEAInd>  &population_front);
void population2front(vector <CCAEAInd>  &mypopulation, vector <vector<double> >  &population_front);
void save_population(vector <CCAEAInd>  &mypopulation, char saveFilename[1024]);

#endif

