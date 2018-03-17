#ifndef __MOEAD_HIGH_GLOBAL_H_
#define __MOEAD_HIGH_GLOBAL_H_

#include "../global.h"
#include "../common.h"
#include "../MOEAD/dmoea.h"
#include "../MOEAD/moeadind.h"
#include "../recomb.h"
#include "../wfg.h"

class TMOEAD_HIGH_GLOBAL
{
public:
	clock_t start,finish;
	clock_t unevolvetime;

	vector <TSOP>  population;  // current population
	CIndividual *indivpoint;    // reference point
	int n_H;
	int  niche;                 // neighborhood size
	int  popsize;                  // population   size
	vector <CIndividual>  ps;
	vector <double> idealpoint;

	TMOEAD_HIGH_GLOBAL(int n_H,int pops,int moead_nei);
	virtual ~TMOEAD_HIGH_GLOBAL();

	void gen_uniformweight(int start_obj_index, int max_value_left, vector<int> coordinate, int sd);
	void init_uniformweight();    // initialize the weights for subproblems

	//分割区域
	vector<vector<int> >lambda;
	//map<int, int> map_index_process;
	vector<vector<int> > sop_partitions;
	vector<set<int> > partition_indexs;
	int maxLayer, maxNiche;

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
	map<int, int> map_weightIndex_popIndex;
	map<int, int> map_popIndex_weightIndex;
	vector<vector<int> > neighbor_overlap_indexs;

	void init_neighborhood();          // calculate the neighborhood of each subproblem
	void init_population();             // initialize the population

	bool update_idealpoint(CIndividual &ind);           // update the approximation of ideal point
	bool update_idealpoint(const double* y_obj);
	void update_problem(CIndividual &child, int id);   // compare and update the neighboring solutions
	void evolution();                                  // mating restriction, recombination, mutation, update

	void execute(int mg, int irun,vector<double>& igd,vector<double> &hv, double &runtime,PTYPE ptype = IDEAL_INDIV);//vector<double>  execute(int sd, int nc, int mg, int rn);          // execute MOEAD

	void save_front(char savefilename[1024]);          // save the pareto front into files
	void save_pos(char saveFilename[1024]);

	void operator=(const TMOEAD_HIGH_GLOBAL &emo);

	void gather_pop_y();
	void gather_populations();
	void calc_qualities(double& igd,double &hv);
	double calc_distance(vector<vector<double> >& ParetoFront);
	double calc_hypervolume(vector<vector<double> >& ParetoFront);

	bool update_idealpoint_from_outside();
	//void sendIdealpoint();
	void sendIdealpoint(double *buf, bool & isSend, int target, MPI_Request* request);
	void update_problem_from_outside();
	void getAndUpdateIndivs(double* b_recv_indiv,int datasize,int source,MPI_Request *request);
	void sendIndiv(double* b_send_indiv,int index,int target,MPI_Request *request,bool &isBegin);
	void update_indiv_from_outside();

	//Parallel Parameter
	//int *pops_sub;
	int *pops_sub_self;
	//int *pops_sub_start;
	double *b_recv,*b_send;
	int s_recv,s_send;

	//更新理想点
	bool isIdealpointUpdate;
	MPI_Request *req_idealpoint_recv;
	MPI_Request *req_idealpoint_send;
	bool *isBegin_Ideal;
	double **recv_idealpoint,**send_idealpoint;
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

#endif
