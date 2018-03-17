#ifndef __MOEAD_GLOBAL_H_
#define __MOEAD_GLOBAL_H_

#include "../global.h"
#include "../common.h"
#include "../MOEAD/moeadind.h"
#include "../recomb.h"
#include "../wfg.h"


class TMOEAD_GLOBAL
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

	TMOEAD_GLOBAL(int n_H,int pops,int moead_nei);
	virtual ~TMOEAD_GLOBAL();

	void gen_uniformweight(int start_obj_index, int max_value_left, vector<int> coordinate, int sd);
	void init_uniformweight();    // initialize the weights for subproblems
	void init_neighborhood();          // calculate the neighborhood of each subproblem
	void init_population();             // initialize the population

	bool update_idealpoint(CIndividual &ind);           // update the approximation of ideal point
	bool update_idealpoint(const double* y_obj, int side);
	void update_problem(CIndividual &child, int id);   // compare and update the neighboring solutions
	void evolution();                                  // mating restriction, recombination, mutation, update

	void execute(int mg, int irun,vector<double>& igd,vector<double> &hv, double &runtime);//vector<double>  execute(int sd, int nc, int mg, int rn);          // execute MOEAD

	void save_front(char savefilename[1024]);          // save the pareto front into files
	void save_pos(char saveFilename[]);

	void operator=(const TMOEAD_GLOBAL &emo);

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
	int *pops_sub;
	int *pops_sub_self;
	int *pops_sub_start;
	double *b_recv,*b_send;
	int s_recv,s_send;

	int prevRank;
	int nextRank;
	//更新理想点
	bool isIdealpointUpdate[2];
	MPI_Request req_idealpoint[4];
	bool isBeginSend[2];
	double *recv_idealpoint[2],*send_idealpoint[2];
	MPI_Status status_idealpoint;

	//更新个体
	bool *recordUpdateofIndiv;
	MPI_Request req_indiv[4];
	bool isLSend;
	bool isRSend;
	int LSize;
	int LLSize;
	int RSize;
	int RRSize;
	double *lrecv_indiv,*lsend_indiv,*rrecv_indiv,*rsend_indiv;
	int s_indiv;
	MPI_Status status_indiv;

	//完成状态
	MPI_Request req_done[4];
	int isDone[4];
	int b_isdone[2];
	MPI_Status status_done;

};

#endif
