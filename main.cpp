/*==========================================================================
//  C++ Implementation of MOEA/D Based on Differential Evolution (DE) for Contest Multiobjective
//  Problems in CEC2009
//
//  Author: Hui Li
//
//  See the details of MOEA/D-DE and test problems in the following papers
//
//  1) H. Li and Q. Zhang, Comparison Between NSGA-II and MOEA/D on a Set of Multiobjective Optimization
//  Problems with Complicated Pareto Sets, Technical Report CES-476, Department of Computer Science,
//  University of Essex, 2007
//
//  2) H. Li and Q. Zhang, Multiobjective Optimization Problems with Complicated Pareto Sets, MOEA/D and NSGA-II,
//  IEEE Transaction on Evolutionary Computation, 2008, to appear.
//
//  If you have any questions about the codes, please contact
//  Dr. Hui Li       at   hzl@cs.nott.ac.uk   or
//  Dr. Qingfu Zhang at   qzhang@essex.ac.uk
//
//  Date: 14/NOV/2008
//
// ===========================================================================*/
//=====================================================*/
//modify by Binson
//Date:July/2016
//=====================================================*/

#include "MOEAD/dmoea.h"
#include "CHEA/cheaalg.h"
#include "CAEA/caeaalg.h"
#include "TestFunction/testInstance.h"
#include "MOEAD - HIGH/dmoea-high.h"
#include "MOEAD - HIGH - GLOBAL/dmoea-high-global.h"
#include "CAEA/ccaeaalg.h"

#include<iostream>
#include<string>
#include<fstream>
#include<vector>
#include<stdlib.h>
using namespace std;

double* hypervolumeRefPoint_TopRight;
double* hypervolumeRefPoint_BottomLeft;

int main(int argc,char* argv[])
{
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);

    printf("the current rank is%d\n",mpi_rank );

//modifiable/>
    char hostname[256];
	int resultLen;
	MPI_Get_processor_name(hostname, &resultLen);

	cout << setprecision(10);
	//sprintf(WORKINGDIR, "%s", "//B8204-4-PC/MPIWorkspace");
	sprintf(WORKINGDIR, "%s", "../");

    cout << "procs : " << numprocs << "    rank : " << mpi_rank <<"    "<<hostname<< endl;

	//int pop;
	char paramFile[256];
	sprintf(paramFile,"%s/labData/param.txt",WORKINGDIR);
	//printf("%s\n",paramFile);
	std::ifstream readf(paramFile);
	//获取相关配置信息
	char tmp[10];
	readf>>tmp;
	if(!strcmp(tmp,"DEBUG")||!strcmp(tmp,"debug"))
		IS_DEBUG = true;
	else if(!strcmp(tmp,"RELEASE") || !strcmp(tmp,"release"))
		IS_DEBUG = false;
	if(IS_DEBUG)
	{
		readf>>tmp;
		if(!strcmp(tmp,"ALL") || !strcmp(tmp,"all"))
		{
			MODE_DEBUG = 1;
			readf>>I_OUT;
		}
		else
			MODE_DEBUG = 0;
	}
	int num_alg=0;
	vector<int> alg;
	readf>>num_alg;
	bool hasMOEAD = false;
	for(int i=0;i<num_alg;i++)
	{
		char tmp[256];
		readf>>tmp;
		if(!strcmp(tmp,"MOEAD_TCH"))
		{
			alg.push_back(0);
			hasMOEAD = true;
		}
		else if(!strcmp(tmp,"MOEAD_PBI"))
		{
			alg.push_back(1);
			hasMOEAD = true;
		}
		else if(!strcmp(tmp,"CAEA"))
			alg.push_back(2);
		else if(!strcmp(tmp,"CHEA"))
			alg.push_back(3);
		else if (!strcmp(tmp, "MOEAD_TCH_HIGH"))
		{
			//cout << "[" << mpi_rank << "] MOEAD_TCH_HIGH" << endl;
			alg.push_back(6);
			hasMOEAD = true;
		}
		else if (!strcmp(tmp, "MOEAD_PBI_HIGH"))
		{
			alg.push_back(7);
			hasMOEAD = true;
		}
		else if (!strcmp(tmp, "MOEAD_TCH_HIGH_GLOBAL"))
		{
			//cout << "[" << mpi_rank << "] MOEAD_TCH_HIGH_GLOBAL" << endl;
			alg.push_back(8);
			hasMOEAD = true;
		}
		else if (!strcmp(tmp, "MOEAD_PBI_HIGH_GLOBAL"))
		{
			alg.push_back(9);
			hasMOEAD = true;
		}
	}
	if(hasMOEAD)
	{
		//MOEAD
		readf>>moead_nei;
	}

	readf >> proc_nei;

//</
    readf>>comm_i_direct;
    if(comm_i_direct>0)
        comm_i_direct_flag=true;
    readf>>comm_i_mig;
    if(comm_i_mig>0)
        comm_i_mig_flag=true;


    ///>

	readf>>total_run;
	int numOfCommunityProblem;
	readf>>numOfCommunityProblem;

	if(mpi_rank == 0)
		cout << "\n[" << mpi_rank << "] -- " << numOfCommunityProblem << " Instances are being tested ---" << endl;;

	char** instances = new char*[numOfCommunityProblem];
	int* numObj = new int[numOfCommunityProblem];
	int* numH = new int[numOfCommunityProblem];
	int* numPop = new int[numOfCommunityProblem];
	int* numGen = new int[numOfCommunityProblem];
    double * px = new double[numOfCommunityProblem];
    double *pm  = new double[numOfCommunityProblem];

    //vector<MGraph>graphs;
    //graphs.resize(numOfCommunityProblem);

	for (int i = 0; i < numOfCommunityProblem; i++)
	{
		instances[i] = new char[256];
		readf >> instances[i];
		readf >> numObj[i];
		readf >> numH[i];
		readf >> numPop[i];
		readf >> numGen[i];
        readf >> px[i];
        readf >> pm[i];
	}
	readf.close();

	for (int inst = 1; inst <= numOfCommunityProblem; inst++)
	{
		// the parameter setting of test instance
        char filename[256];
	    sprintf(filename,"%s/graph model/%s",WORKINGDIR,instances[inst-1]);

        MGraph *G = new MGraph;
        init_graphModel(G,filename);
		nvar = G->numVertexes;
		nobj = numObj[inst - 1];
		n_H =  numH[inst-1];
		pops = numPop[inst - 1];
        sprintf(strTestInstance, "%s", instances[inst - 1]);
		total_gen = numGen[inst - 1];
        graph_xcross_p=px[inst-1];
        graph_mutation_p=pm[inst-1];


//</modifiable
    //printf("----usage:%s --is_Debug(Yes,option:Re) --paraType(IDEAL_INDIV) --paraTopology(FULL)\n", argv[0]);
    for (int i = 1; i < argc; i++) {
        if(strcmp(argv[i],"-II")==0)paraType=IDEAL_INDIV;
        if(strcmp(argv[i],"-ID")==0)paraType=IDEAL;
        if(strcmp(argv[i],"-IV")==0)paraType=INDIV;
        if(strcmp(argv[i],"-F")==0)paraTopology=FULL;
        if(strcmp(argv[i],"-B")==0)paraTopology=BI_RING;
        //if(strcmp(argv[i],"-S")==0)speedups=true;
    }
//modifiable/>
		//printf("\n -- Instance: %s, %d variables %d objectives \n\n", strTestInstance, nvar, nobj);
		int algno=0;
		for(int alg_index = 0; alg_index < alg.size(); alg_index++)
		{
			algno = alg[alg_index];

			std::fstream fouttime;
            std::fstream foutqmetric;
			//std::fstream fouthv;
			//char hvfile[1024];
            char   qfile[1024];
			char   timefile[1024];
			if(mpi_rank == 0)
			{//只在主机输出
				if(IS_DEBUG)
				{
					//
					switch (algno)
					{
						case 0:
						{
							sprintf(timefile, "%s/labData/LOG/LOG_MOEAD_TCH_%s(%d)_%d_%d", WORKINGDIR, strTestInstance, nobj, pops, total_gen);
							sprintf(qfile, "%s/labData/Q/Q_MOEAD_TCH_%s(%d)_%d_%d", WORKINGDIR, strTestInstance, nobj, pops, total_gen);
							//sprintf(hvfile, "%s/labData/HV/HV_MOEAD_TCH_%s(%d)_%d_%d", WORKINGDIR, strTestInstance, nobj, pops, total_gen);
							cout << "[" << mpi_rank << "] MOEA/D-TCH--Inst: " << strTestInstance << " (nvar: " << nvar << " , nobj: " << nobj << " )--pops: " << pops << " --total_gen:" << total_gen << endl;
							break;
						}
						case 1:
						{
							sprintf(timefile, "%s/labData/LOG/LOG_MOEAD_PBI_%s(%d)_%d_%d", WORKINGDIR, strTestInstance, nobj, pops, total_gen);
							sprintf(qfile, "%s/labData/Q/Q_MOEAD_PBI_%s(%d)_%d_%d", WORKINGDIR, strTestInstance, nobj, pops, total_gen);
							//sprintf(hvfile, "%s/labData/HV/HV_MOEAD_PBI_%s(%d)_%d_%d", WORKINGDIR, strTestInstance, nobj, pops, total_gen);
							cout << "[" << mpi_rank << "] MOEA/D-PBI --Inst: " << strTestInstance << " (nvar: " << nvar << " , nobj: " << nobj << " )--pops: " << pops << " --total_gen:" << total_gen << endl;
							break;
						}
						case 2:
						{
							sprintf(timefile, "%s/labData/LOG/LOG_CAEA_%s(%d)_%d_%d", WORKINGDIR, strTestInstance, nobj, pops, total_gen);
							sprintf(qfile, "%s/labData/Q/Q_CAEA_%s(%d)_%d_%d", WORKINGDIR, strTestInstance, nobj, pops, total_gen);
							//sprintf(hvfile, "%s/labData/HV/HV_CAEA_%s(%d)_%d_%d", WORKINGDIR, strTestInstance, nobj, pops, total_gen);
							cout << "[" << mpi_rank << "] CAEA--Inst: " << strTestInstance << " (nvar: " << nvar << " , nobj: " << nobj << " )--pops: " << pops << " --total_gen:" << total_gen << endl;
							break;
						}
						case 3:
						{
							sprintf(timefile, "%s/labData/LOG/LOG_CHEA_%s(%d)_%d_%d", WORKINGDIR, strTestInstance, nobj, pops, total_gen);
							sprintf(qfile, "%s/labData/Q/Q_CHEA_%s(%d)_%d_%d", WORKINGDIR, strTestInstance, nobj, pops, total_gen);
							//sprintf(hvfile, "%s/labData/HV/HV_CHEA_%s(%d)_%d_%d", WORKINGDIR, strTestInstance, nobj, pops, total_gen);
							cout << "[" << mpi_rank << "] CHEA--Inst: " << strTestInstance << " (nvar: " << nvar << " , nobj: " << nobj << " )--pops: " << pops << " --total_gen:" << total_gen << endl;
							break;
						}
						case 4:
						{
							sprintf(timefile, "%s/labData/LOG/LOG_MOEAD_TCH_GLOBAL_%s(%d)_%d_%d", WORKINGDIR, strTestInstance, nobj, pops, total_gen);
							sprintf(qfile, "%s/labData/Q/Q_MOEAD_TCH_GLOBAL_%s(%d)_%d_%d", WORKINGDIR, strTestInstance, nobj, pops, total_gen);
							//sprintf(hvfile, "%s/labData/HV/HV_MOEAD_TCH__GLOBAL%s(%d)_%d_%d", WORKINGDIR, strTestInstance, nobj, pops, total_gen);
							cout << "[" << mpi_rank << "] MOEA/D-TCH-GLOBAL--Inst: " << strTestInstance << " (nvar: " << nvar << " , nobj: " << nobj << " )--pops: " << pops << " --total_gen:" << total_gen << endl;
							break;
						}
						case 5:
						{
							sprintf(timefile, "%s/labData/LOG/LOG_MOEAD_PBI_GLOBAL_%s(%d)_%d_%d", WORKINGDIR, strTestInstance, nobj, pops, total_gen);
							sprintf(qfile, "%s/labData/Q/Q_MOEAD_PBI_GLOBAL%s(%d)_%d_%d", WORKINGDIR, strTestInstance, nobj, pops, total_gen);
							//sprintf(hvfile, "%s/labData/HV/HV_MOEAD_PBI_GLOBAL%s(%d)_%d_%d", WORKINGDIR, strTestInstance, nobj, pops, total_gen);
							cout << "[" << mpi_rank << "] MOEA/D-PBI-GLOBAL--Inst: " << strTestInstance << " (nvar: " << nvar << " , nobj: " << nobj << " )--pops: " << pops << " --total_gen:" << total_gen << endl;
							break;
						}
						case 6:
						{
							sprintf(timefile, "%s/labData/LOG/LOG_MOEAD_TCH_HIGH_%s(%d)_%d_%d", WORKINGDIR, strTestInstance, nobj, pops, total_gen);
							sprintf(qfile, "%s/labData/Q/Q_MOEAD_TCH_HIGH_%s(%d)_%d_%d", WORKINGDIR, strTestInstance, nobj, pops, total_gen);
							//sprintf(hvfile, "%s/labData/HV/HV_MOEAD_TCH_HIGH_%s(%d)_%d_%d", WORKINGDIR, strTestInstance, nobj, pops, total_gen);
							cout << "[" << mpi_rank << "] MOEA/D-TCH_HIGH--Inst: " << strTestInstance << " (nvar: " << nvar << " , nobj: " << nobj << " )--pops: " << pops << " --total_gen:" << total_gen << endl;
							break;
						}
						case 7:
						{
							sprintf(timefile, "%s/labData/LOG/LOG_MOEAD_PBI_HIGH_%s(%d)_%d_%d", WORKINGDIR, strTestInstance, nobj, pops, total_gen);
							sprintf(qfile, "%s/labData/Q/Q_MOEAD_PBI_HIGH_%s(%d)_%d_%d", WORKINGDIR, strTestInstance, nobj, pops, total_gen);
							//sprintf(hvfile, "%s/labData/HV/HV_MOEAD_PBI_HIGH_%s(%d)_%d_%d", WORKINGDIR, strTestInstance, nobj, pops, total_gen);
							cout << "[" << mpi_rank << "] MOEA/D-PBI_HIGH --Inst: " << strTestInstance << " (nvar: " << nvar << " , nobj: " << nobj << " )--pops: " << pops << " --total_gen:" << total_gen << endl;
							break;
						}
						case 8:
						{
							sprintf(timefile, "%s/labData/LOG/LOG_MOEAD_TCH_HIGH_GLOBAL_%s(%d)_%d_%d", WORKINGDIR, strTestInstance, nobj, pops, total_gen);
							sprintf(qfile, "%s/labData/Q/Q_MOEAD_TCH_HIGH_GLOBAL_%s(%d)_%d_%d", WORKINGDIR, strTestInstance, nobj, pops, total_gen);
							//sprintf(hvfile, "%s/labData/HV/HV_MOEAD_TCH_HIGH_GLOBAL_%s(%d)_%d_%d", WORKINGDIR, strTestInstance, nobj, pops, total_gen);
							cout << "[" << mpi_rank << "] MOEA/D-TCH_HIGH_GLOBAL--Inst: " << strTestInstance << " (nvar: " << nvar << " , nobj: " << nobj << " )--pops: " << pops << " --total_gen:" << total_gen << endl;
							break;
						}
						case 9:
						{
							sprintf(timefile, "%s/labData/LOG/LOG_MOEAD_PBI_HIGH_GLOBAL_%s(%d)_%d_%d", WORKINGDIR, strTestInstance, nobj, pops, total_gen);
							sprintf(qfile, "%s/labData/Q/Q_MOEAD_PBI_HIGH_GLOBAL_%s(%d)_%d_%d", WORKINGDIR, strTestInstance, nobj, pops, total_gen);
							//sprintf(hvfile, "%s/labData/HV/HV_MOEAD_PBI_HIGH_GLOBAL_%s(%d)_%d_%d", WORKINGDIR, strTestInstance, nobj, pops, total_gen);
							cout << "[" << mpi_rank << "] MOEA/D-PBI_HIGH_GLOBAL --Inst: " << strTestInstance << " (nvar: " << nvar << " , nobj: " << nobj << " )--pops: " << pops << " --total_gen:" << total_gen << endl;
							break;
						}
					}
                    sprintf(timefile, "%s_%dnp_%s_%s",timefile,numprocs,getElemTypeName(paraType),getElemTypeName(paraTopology));
                    sprintf(qfile, "%s_%dnp_%s_%s",qfile,numprocs,getElemTypeName(paraType),getElemTypeName(paraTopology));
                    //sprintf(hvfile, "%s_%dnp_%s_%s",hvfile,numprocs,getElemTypeName(paraType),getElemTypeName(paraTopology));



                    sprintf(timefile, "%s.dat",timefile);
					sprintf(qfile, "%s.dat",qfile);
					//sprintf(hvfile, "%s.dat",hvfile);

                    fouttime.open(timefile, std::ios::out);
					fouttime << setprecision(10);
					foutqmetric.open(qfile, std::ios::out);
					foutqmetric << setprecision(10);
					//fouthv.open(hvfile, std::ios::out);
					//fouthv << setprecision(10);
				}

			}
            cout<<"px: "<<graph_xcross_p<<" pm: "<<graph_mutation_p<<endl;
			for (int run = 1; run <= total_run; run++)
			{
				if(mpi_rank == 0)
					cout<<"-- "<<run<<"-th run  -- "<<endl;
				//同步
				//cout<<"["<<mpi_rank<<"]"<<run<<"-"<<nvar<<"-"<<nobj<<endl;

				MPI_Barrier(MPI_COMM_WORLD);
				//vector<double> hv;
				vector<double> qmetric;
				double runtime = 0;
				seed = int(clock()) % 1377;
				rnd_uni_init = -(long)seed;

				switch(algno)
				{
					//case 0:	//MOEAD_TCH
						//{
							//tDecompType = _Decomposition_TCH1;
							//sprintf(STR_MOEAD,"MOEAD_TCH");
							//TMOEAD MOEAD(n_H,pops, moead_nei);
							//MOEAD.execute(total_gen, run, qmetric, runtime, paraType);
						//}
						//break;
					//case 1:	//MOEAD_PBI
						//{
							//tDecompType = _Decomposition_PBI;
							//sprintf(STR_MOEAD,"MOEAD_PBI");
							//TMOEAD MOEAD(n_H,pops, moead_nei);
							//MOEAD.execute(total_gen, run, qmetric, runtime, paraType);

						//}
						//break;
					case 2: //CAEA
						{
							//CCAEA caea(pops);
							//caea.execute(total_gen, run, qmetric, runtime, paraType);
							CAEA caea(pops, n_H);
							caea.execute(total_gen, run, qmetric, runtime, paraType,G);

						}
						break;
					case 3: //CHEA
						{
							CHEA chea(pops, n_H);
							chea.execute(total_gen, run, qmetric, runtime,paraType,G);
						}
						break;
						/*
					case 4://MOEAD TCH GLOBAL
						{
							tDecompType = _Decomposition_TCH1;
							sprintf(STR_MOEAD, "MOEAD_TCH");
							TMOEAD_GLOBAL MOEAD(n_H, pops, moead_nei);
							MOEAD.execute(total_gen, run, qmetric, runtime);

						}
						break;
					case 5:	//MOEAD_PBI	 GLOBAL
						{
							tDecompType = _Decomposition_PBI;
							sprintf(STR_MOEAD, "MOEAD_PBI");
							TMOEAD_GLOBAL MOEAD(n_H, pops, moead_nei);
							MOEAD.execute(total_gen, run, qmetric, runtime);

						}
						break;
						*/
					//case 6:	//MOEAD_TCH_HIGH
					//{
						////cout << "[" << mpi_rank << "] 1" << endl;
						//tDecompType = _Decomposition_TCH1;

						//sprintf(STR_MOEAD, "MOEAD_TCH_HIGH");
						////cout << "[" << mpi_rank << "]" << STR_MOEAD << endl;
						//TMOEAD_HIGH MOEAD(n_H, pops, moead_nei);
						//MOEAD.execute(total_gen, run, qmetric, runtime, paraType);

					//}
					//break;
					//case 7:	//MOEAD_PBI_HIGH
					//{
						//tDecompType = _Decomposition_PBI;
						//sprintf(STR_MOEAD, "MOEAD_PBI_HIGH");
						//TMOEAD_HIGH MOEAD(n_H, pops, moead_nei);
						//MOEAD.execute(total_gen, run, qmetric, runtime, paraType);

					//}
					//break;
					//case 8:	//MOEAD_TCH_HIGH_GLOBAL
					//{
						////cout << "[" << mpi_rank << "] 1" << endl;
						//tDecompType = _Decomposition_TCH1;

						//sprintf(STR_MOEAD, "MOEAD_TCH_HIGH_GLOBAL");
						////cout << "[" << mpi_rank << "]" << STR_MOEAD << endl;
						//TMOEAD_HIGH_GLOBAL MOEAD(n_H, pops, moead_nei);
						//MOEAD.execute(total_gen, run, qmetric, runtime, paraType);

					//}
					//break;
					//case 9:	//MOEAD_PBI_HIGH_GLOBAL
					//{
						//tDecompType = _Decomposition_PBI;
						//sprintf(STR_MOEAD, "MOEAD_PBI_HIGH_GLOBAL");
						//TMOEAD_HIGH_GLOBAL MOEAD(n_H, pops, moead_nei);
						//MOEAD.execute(total_gen, run, qmetric, runtime, paraType);
					//}
					//break;
				}

				if(mpi_rank ==0)
				{//保存Q和HV数据
					//cout << "==========="<<run<<"=============" << endl;
					for(int k = 0; k < qmetric.size(); k++)
					{
						if (k%2 == 0)
							foutqmetric << qmetric[k] << " ";
						else
						{
							foutqmetric << qmetric[k] << '\n';
						}
					}
					//foutqmetric << "\n";
					qmetric.clear();

					//for(int k =0;k < hv.size();k++)
					//{
						//if(k%2 == 0)
							//fouthv<<hv[k]<<" ";
						//else
							//fouthv<<hv[k]<<'\n';
					//}
					//fouthv<<'\n';
					//hv.clear();

					fouttime << run << " " << runtime << "\n";
				}
			}
			//fouthv.close();
			fouttime.close();
			foutqmetric.close();
		}
		hv_ref_lb.clear();
		hv_ref_ub.clear();
	}
	delete []instances;
	delete []numObj;
	delete []numGen;
	delete []numPop;
	alg.clear();
	MPI_Finalize();
}
