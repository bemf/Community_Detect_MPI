#ifndef __GLOBAL_H_
#define __GLOBAL_H_

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <memory.h>
#include <vector>
#include <map>
#include <cassert>
#include <algorithm>
#include "mpi.h"
#include "random.h"
#include <iomanip>
#include <set>
#include <queue>

using namespace std;

//#include "cec09.h"

#define DEBUG 1

#define PI  3.1415926535897932384626433832795

#define MAXVEX 5000
#define INFINITY 65535
typedef char VertexType;
typedef int EdgeType;
typedef struct
{
	VertexType vexs[MAXVEX];
	EdgeType arc[MAXVEX][MAXVEX];

    int Adjacency[MAXVEX][MAXVEX];
    int degree[MAXVEX];
	int numVertexes, numEdges;
}MGraph;
void init_graphModel(MGraph*,char*);
//------------- Parameters in test instance ------------------

extern int     nvar,  nobj;                    //  the number of variables and objectives
extern int n_H;

//double  lowBound = 0,   uppBound = 1;   //  lower and upper bounds of variables


//void obj_eval_ywq(vector <double>& x_var, vector <double>& y_obj);

enum TDecompositionType {_Decomposition_TCH1, _Decomposition_TCH2, _Decomposition_PBI};

extern TDecompositionType    tDecompType ;
extern char STR_MOEAD[256];

/** Possible relations between points in multiple objective space */
enum TCompare {_Pareto_Dominating, _Pareto_Dominated, _Pareto_Nondominated, _Pareto_Equal};


extern void (*pFunc[])(double *, double *, const unsigned int);

extern vector<double>lowBound,uppBound;

extern vector<double> hv_ref_ub;
extern vector<double> hv_ref_lb;
extern double cubeHV;

extern char    strTestInstance[256];


//******** Common parameters in MOEAs **********************************************
extern int		total_gen,    //  the maximal number of generations
		total_run,//1,      //  the maximal number of runs
		pops ;    //  the population size

extern  int I_OUT ;
extern int I_UPDATE;
	//	nfes;             //  the number of function evluations
//**********************************************************************************

extern int moead_nei ;
extern bool IS_DEBUG;
extern int MODE_DEBUG;
//------------- Parameters in random number ------------------
extern int     seed ;
extern long    rnd_uni_init;

extern int		etax, etam;   // distribution indexes of crossover and mutation

//extern double  realx,  realm,  realb;    // crossover, mutation, selection probabilities

/******************************************************
* Parameters for Parallel
*/
extern int comm_i;
extern char WORKINGDIR[1024];
extern int proc_nei;
extern int comm_i_direct,comm_i_mig;
extern bool comm_i_direct_flag,comm_i_mig_flag;

extern int numprocs;
extern int mpi_rank;
extern int root_rank;
/*******************************************************/
extern const int TAG_IDEALPOINT;
extern const int TAG_ANCHORPOINT ;
extern const int TAG_INDIV ;
extern const int TAG_DONE ;
extern const int RECV ;
extern const int SEND;
extern const int L_RECV;
extern const int L_SEND;
extern const int R_RECV;
extern const int R_SEND;
extern const int LEFT;
extern const int RIGHT;


extern bool is_first;

enum PTYPE {IDEAL,INDIV,IDEAL_INDIV,NONE};
extern PTYPE paraType;

inline char *getElemTypeName(PTYPE type)
{
    switch (type) {
        case IDEAL: return "IDEAL";
        case INDIV: return "INDIV";
        case IDEAL_INDIV: return "IDEAL_INDIV";
        case NONE: return "NONE";
        default: return "error";
    }
}

enum TOPOLOGY {DYNAMIC,UNI_RING,BI_RING,FULL,NEI};
extern TOPOLOGY paraTopology;
inline char *getElemTypeName(TOPOLOGY type)
{
    switch (type) {
        case UNI_RING: return "UNI_RING";
        case BI_RING: return "BI_RING";
        case FULL: return "FULL";
        case NEI: return "NEI";
        case DYNAMIC: return "DYNAMIC";
        default: return "error";
    }
}

inline void mpi_error_handler(int error_code,int my_rank,const char* file,const char* function,int line)
{
   if (error_code != MPI_SUCCESS)
   {

       char error_string[1024];
       int length_of_error_string, error_class;

       MPI_Error_class(error_code, &error_class);
       MPI_Error_string(error_class, error_string, &length_of_error_string);
       fprintf(stderr,"%s(%d)-<%s>:",file,line,function);
       fprintf(stderr, "%3d: %s\n", my_rank, error_string);
       MPI_Error_string(error_code, error_string, &length_of_error_string);
       fprintf(stderr,"%s(%d)-<%s>:",file,line,function);
       fprintf(stderr, "%3d: %s\n", my_rank, error_string);
       MPI_Abort(MPI_COMM_WORLD, error_code);
    }
}
extern bool speedups;
inline void idle()
{
    long int times=1000000;
    long int i=0;
    for(;i<times;i++)i=i;
}
extern int err_code;

extern double graph_xcross_p, graph_mutation_p;
#endif
