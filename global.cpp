#include "global.h"

//------------- Parameters in test instance ------------------

int     nvar,  nobj;                    //  the number of variables and objectives
int n_H;

//double  lowBound = 0,  uppBound = 1;  //  lower and upper bounds of variables
vector<double>  lowBound,uppBound;

vector<double> hv_ref_ub;
vector<double> hv_ref_lb;
double cubeHV;

char    strTestInstance[256];

TDecompositionType    tDecompType =_Decomposition_PBI;////  _Decomposition_TCH1;//
char STR_MOEAD[256];
//******** Common parameters in MOEAs **********************************************
int		total_gen = 100,//250,//250, //500 1500 1000,//500,//    //  the maximal number of generations
		total_run = 10, //30,//      //  the maximal number of runs
		pops   = 100;//200,    //  the population size
 int I_OUT = 25;
 int I_UPDATE = 10;
		//nfes;             //  the number of function evluations
//**********************************************************************************

int moead_nei = 20;
bool IS_DEBUG = true;
int MODE_DEBUG = 0;
//------------- Parameters in random number ------------------
int     seed    = 177;
long    rnd_uni_init;

///15 20
int		etax = 20, etam = 20;   // distribution indexes of crossover and mutation

/******************************************************
* Parameters for Parallel
*/
char WORKINGDIR[1024];
int proc_nei = 10;

int comm_i=10;
int comm_i_direct = 10;
int comm_i_mig=10;
bool comm_i_direct_flag=false;
bool comm_i_mig_flag=false;

int numprocs;
int mpi_rank;
int root_rank = 0;
/*******************************************************/

const int TAG_IDEALPOINT = 1001;
const int TAG_ANCHORPOINT = 1002;
const int TAG_INDIV = 2000;
const int TAG_DONE = 9999;
const int RECV = 0;
const int SEND = 1;
const int L_RECV = 0;
const int L_SEND = 1;
const int R_RECV = 2;
const int R_SEND = 3;
const int LEFT = 0;
const int RIGHT = 1;

PTYPE paraType =IDEAL_INDIV;// IDEAL;// INDIV; //
TOPOLOGY paraTopology =   BI_RING;//
//TOPOLOGY paraTopology =  FULL;//NEI;// BI_RING;//
bool is_first=true;
int err_code;

void init_graphModel(MGraph *G,char*filename)
{
	std::ifstream infile(filename);
	std::ifstream infile1(filename);
	string temp;
	char *s;
	int numv = 0, nume = 0;
	while (getline(infile1, temp))
	{
		numv++;
	}
	cout << "filename: "<< filename << endl;
	cout << "numV: " << numv << endl;
	G->numVertexes = numv;
	for (int i = 0; i < numv; ++i)
	{
		for (int j = 0; j < numv; ++j){
			G->arc[i][j] = 0;
		}
           G->degree[i]=0;
	}
	bool isvertex = false;
	vector <int>  tempvertex;
	vector <int> tempedges;
    int numEdges = 0;
	while (getline(infile, temp))
	{
		isvertex = true;
		s = (char *)temp.data();
		const char *d = " ";
		char *p;
		p = strtok(s, d);
		while (p){
			int b = atoi(p);//将char*类型转换为int
			if (isvertex){
				tempvertex.push_back(b);
				isvertex = false;
			}
			else{
				tempedges.push_back(b);
			}
			//cout << p << endl;
			p = strtok(NULL, d);
		}
		int i, j;
		for (vector<int>::iterator it = tempedges.begin(); it != tempedges.end();)
		{
			i = tempvertex.size() - 1;
			j = *(it);
			G->arc[i][j] = 1;
            numEdges++;
            G->Adjacency[i][G->degree[i]++]=j;
			it = tempedges.erase(it);
		}
	}
    G->numEdges=numEdges/2;
}



double graph_xcross_p, graph_mutation_p;
