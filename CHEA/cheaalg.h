#ifndef __CHEA_CLASS_H_
#define __CHEA_CLASS_H_

#include "../global.h"
#include "../common.h"
#include "../recomb.h"
#include "../wfg.h"
#include "cheaind.h"

class CHEA
{
public:
	CHEA(int pop_size, int hype_intercept);
	~CHEA();
	//void execute(int mg, int irun, vector<double>& igd, vector<double>& hv, double &runtime, PTYPE ptype);

	void evolution_tour_select();
	void evolution_tour_select_2();
	vector<int> subp_index_on_edge;
	int size_subp_on_edge;

	void init_uniformweight();
	void gen_uniformweight(int start_obj_index, int max_value_left, vector<int> coordinate, int H);
	void init_population();
	//基于社区检测的锦标赛算法选择
	void evolution_tour_select_2(MGraph *G);
	//基于社区检测的初始化
	void init_population(MGraph *G);
	//社区检测执行
	void execute(int mg, int irun, vector<double>& qmetric, double &runtime, PTYPE ptype, MGraph *G);

	void gen_SubProblems(int start_obj_index, int max_value_left, vector<int> coordinate, int &index_count);

	inline void init_neighborhood()//为每一个子问题计算邻居子问题编号
	{
		for (int i = 0; i < popsize; i++)
		{
			//population[i].get_vicinity(hyperplane_Intercept);
			population[i].neighborhood = gen_vicinities(population[i].V_obj_o ,hyperplane_Intercept);
			/*cout << "[" << mpi_rank << "] (" << i << ") ";
			for (int j = 0; j < population[i].neighborhood.size(); j++)
				cout << population[i].neighborhood[j] << " ";
			cout << endl;
			*/
		}
	}

	void operator=(const CHEA &emo);

	void gather_pop_y();
	void gather_populations();
	void calc_qualities(double& igd, double &hv);
    void calc_qmetric(double&,int&,MGraph*);
	double calc_distance(vector<vector<double> >& ParetoFront);
	double calc_hypervolume(vector<vector<double> >& ParetoFront);

	vector <CHEA_SOP>  population;
	vector <CHEAInd>  ps;
	CHEAInd onechild;

	int popsize;
	int hyperplane_Intercept;
	double per_intercept;

	vector<vector <double> > AnchorPoint;
	vector <double> TrueNadirPoint;
	vector <double> IdealPoint;
	vector <double> ReferencePoint;

	virtual bool Pareto_HyperVolume_compare_sectorialgrid(const CHEAInd& ind ,CHEAInd& replacedInd);

	bool update_extreme_point(CHEAInd& ind);
	bool update_ideal_point(const double* y_obj);

	void update_partition();

	virtual double GetHyperVolume(const CHEAInd&  ind, vector <double> &ref_ponit_cal);

	int  tour_selection_hv(vector <CHEA_SOP>  &population);
	int  tour_selection_Q(vector <CHEA_SOP>  &population,MGraph*G);
	double calc_q(int index,vector <CHEA_SOP>  &population,MGraph*G);
	double  tour_selection_hv_difference(int p, vector <CHEA_SOP>  &population);

	void save_front(char savefilename[1024]);          // save the pareto front into files
	void save_pos(char saveFilename[1024]);

	clock_t start, finish;
	clock_t unevolvetime;

	//分割区域
	vector<vector<int> >lambda;
	//map<int, int> map_index_process;
	vector<vector<int> > sop_partitions;
	vector<set<int> > partition_indexs;
	int maxLayer;

	bool lambda_compare(vector<int> &l1, vector<int> & l2);
	void allocateForProcess();
	void allocateForProcess_Linear();
	vector<set<int> > partition_overlap_indexs;
	void partition_overlap(int maxLayer);

	set<int> neighbor_partitions_set;
	vector<int>neighbor_partitions;
	void partition_neighbor();
	int partition_niche;
	int overlap_start_index;

	bool update_idealpoint_from_outside();
	void sendIdealpoint(double *buf, bool &isSend, int target, MPI_Request* request);
	void getAndUpdateIndivs(double* b_recv_indiv, int datasize, int source, MPI_Request *request);
	void update_indiv_from_outside();


	//Parallel Parameter
	//int *pops_sub;
	int *pops_sub_self;
	//int *pops_sub_start;
	double *b_recv, *b_send;
	int s_recv, s_send;

    double *c_recv,*c_send;
    int cs_recv,cs_send;


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
    double **indiv_buf_direct;
	int s_indiv;
	int indivBufSize;
	MPI_Status status_indiv;

    vector<CHEAInd> send_indivs;

	//完成状态
	MPI_Request *req_done_recv;
	MPI_Request *req_done_send;
	int *isDoneRecv;
	int *isDoneSend;
	int *b_isdone;
	MPI_Status status_done;
};

void population2front(vector <CHEA_SOP>  &mypopulation, vector <CHEAInd>  &population_front);
void population2front(vector <CHEA_SOP>  &mypopulation, vector <vector<double> >  &population_front);
void save_population(vector <CHEAInd>  &mypopulation, char saveFilename[1024]);



#endif
