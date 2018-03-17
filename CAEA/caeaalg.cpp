

#ifndef __CEC09_H_
//#include "cec09.h"
#endif

#include "caeaalg.h"

CCAEA::CCAEA(int pop_size){
	start = clock();
	popsize = pop_size;

	//计算局部进化区域分配大小
	pops_sub_self = new int[numprocs];
	pops_sub_start = new int[numprocs];
	int mod_tmp = popsize %numprocs;
	int v_tmp = popsize / numprocs;
	if (mod_tmp == 0)
		pops_sub_self[0] = v_tmp;
	else
		pops_sub_self[0] = v_tmp + 1;
	pops_sub_start[0] = 0;

	for (int i = 1; i<numprocs; i++)
	{
		pops_sub_self[i] = v_tmp + (mod_tmp > i ? 1 : 0);
		pops_sub_start[i] = pops_sub_start[i - 1] + pops_sub_self[i - 1];
	}

	int tmp = nvar + nobj;
	s_send = tmp * pops_sub_self[mpi_rank];
	s_recv = tmp * popsize*numprocs;
	b_send = new double[s_send];
	if (mpi_rank == 0)
		b_recv = new double[s_recv];
	cout << "[" << mpi_rank << "] start:" << pops_sub_start[mpi_rank] << "  size: " << pops_sub_self[mpi_rank]<<"send_size: "<<s_send<<"recv_size : "<<s_recv << endl;


	rcounts_y = new int[numprocs];
	for (int i = 0; i<numprocs; i++)
		rcounts_y[i] = nobj *pops_sub_self[i];
	displs_y = new int[numprocs];
	displs_y[0] = 0;
	for (int i = 1; i<numprocs; i++)
		displs_y[i] = displs_y[i - 1] + rcounts_y[i - 1];

	rcounts_a = new int[numprocs];
	for (int i = 0; i<numprocs; i++)
		rcounts_a[i] = (nvar + nobj) *pops_sub_self[i];
	displs_a = new int[numprocs];
	displs_a[0] = 0;
	for (int i = 1; i<numprocs; i++)
		displs_a[i] = displs_a[i - 1] + rcounts_a[i - 1];

}

CCAEA::~CCAEA(){
	delete b_send;
	if (mpi_rank == 0)
		delete b_recv;
	population.clear();
	ps.clear();
	TrueNadirPoint.clear();
	IdealPoint.clear();
	ReferencePoint.clear();
	delete rcounts_y;
	rcounts_y = NULL;
	delete displs_y;
	displs_y = NULL;
	delete rcounts_a;
	rcounts_a = NULL;
	delete displs_a;
	displs_a= NULL;
}


bool CCAEA::update_extreme_point(CCAEAInd& ind)
{	 /*
	bool bAnchorUpdated = false;
	bool bTrueNadirUpdated = false;

	int j;

	for( j=0;j<nobj ;j++)
	{
		if( ind.y_obj[j] < AnchorPoint[j][j] || (ind.y_obj[j] == AnchorPoint[j][j] && ind.y_obj[(j+1)%2] < AnchorPoint[j][(j+1)%2] )  )
		//if(ind.y_obj[j] < ObserverPoint[j] )
		{
			bAnchorUpdated = true;
			//ObserverPoint[j] = ind.y_obj[j];
			AnchorPoint[j][j] = ind.y_obj[j];
			AnchorPoint[j][(j+1)%2] = ind.y_obj[(j+1)%2];
		}
	}

	for( j=0;j<nobj ;j++)
	{
		if( ind.y_obj[j] >TrueNadirPoint[j] )
		{
			bTrueNadirUpdated = true;
			TrueNadirPoint[j] = ind.y_obj[j];
		}
	}

	//if(bAnchorUpdated || bTrueNadirUpdated)
	//{
	//	for( j=0;j<nobj ;j++)
	//	{
	//			ReferencePoint[j] = ( TrueNadirPoint[j] + 1e3 * (TrueNadirPoint[j] - ObserverPoint[j]) );
	//	}
	//}

	return bAnchorUpdated;
	*/
	bool bIdealUpdated = false;
	bool bTrueNadirUpdated = false;
	bool bUpdate_item = false;
	int j;

	for (j = 0; j < nobj; j++)
	{
		if (ind.y_obj[j] < IdealPoint[j])
		{
			bIdealUpdated = true;
			IdealPoint[j] = ind.y_obj[j];
			bUpdate_item = true;
		}
		if (ind.y_obj[j] > TrueNadirPoint[j])
		{
			bTrueNadirUpdated = true;
			TrueNadirPoint[j] = ind.y_obj[j];
		}
		if (bUpdate_item || bTrueNadirUpdated)
		{
			ReferencePoint[j] = (TrueNadirPoint[j]) + 1e3 *(TrueNadirPoint[j] - IdealPoint[j]);
			bUpdate_item = false;
			bTrueNadirUpdated = false;
		}
	}

	isIdealpointUpdate |= bIdealUpdated;

	return bIdealUpdated;
}

bool CCAEA::update_ideal_point(const double* y_obj )
{
	bool bIdealUpdated = false;
	bool bUpdate_item = false;
	int j;

	for (j = 0; j < nobj; j++)
	{
		if (y_obj[j] < IdealPoint[j])
		{
			bIdealUpdated = true;
			IdealPoint[j] = y_obj[j];
			bUpdate_item = true;
		}

		if (bUpdate_item)
		{
			ReferencePoint[j] = (TrueNadirPoint[j]) + 1e3 *(TrueNadirPoint[j] - IdealPoint[j]);
			bUpdate_item = false;
		}
	}

	isIdealpointUpdate |= bIdealUpdated;

	return bIdealUpdated;
}

void CCAEA::init_population()
{
	vector <bool>  sectorpopFlag;

	vector <CCAEAInd>  initial_indi_pop;
	CCAEAInd rest_ind;
	CCAEAInd ind;

	for(int n=0;n<popsize;n++)//统计极值点or观测点  初始个体数目设置
	{
		sectorpopFlag.push_back(false);
		ind.rnd_init();
		//nfes++;
		initial_indi_pop.push_back(ind);//insert_sectorialgrid(ind, n);
	}

	//initial_observation_and_reference_point(initial_indi_pop[0]);

	for(int j=0;j<nobj ;j++)
	{
		IdealPoint.push_back(initial_indi_pop[0].y_obj[j]);
		TrueNadirPoint.push_back(initial_indi_pop[0].y_obj[j]);
		ReferencePoint.push_back(initial_indi_pop[0].y_obj[j]);
		//AnchorPoint.push_back( initial_indi_pop[0].y_obj );//ObserverPoint.push_back( initial_indi_pop[0].y_obj[j] );//
		//TrueNadirPoint.push_back( initial_indi_pop[0].y_obj[j] );//ReferencePoint.push_back( initial_indi_pop[0].y_obj[j] );
	}
	//for(int j=0;j<nobj ;j++)
	//{
	//	ObserverPoint.push_back( AnchorPoint[j][j] );
	//	//PseudoNadirPoint.push_back( AnchorPoint[(j+1)%2][j] );
	//}
	//for(int j=0;j<nobj ;j++)
	//{
	//	ReferencePoint.push_back( TrueNadirPoint[j] + 1e3 * (TrueNadirPoint[j] - ObserverPoint[j]) ) ;
	//}

	for(int n=1;n<popsize;n++)//统计极值点or观测点  初始个体数目设置
	{
		update_extreme_point(initial_indi_pop[n]);
	}

	for(int n=0;n<popsize;n++)//统计极值点or观测点  初始个体数目设置
	{
		population.push_back(initial_indi_pop[0]);
	}

	vector <CCAEAInd>  initial_rest_indi_pop;
	for(int n=0;n< popsize;n++)//计算角度 第一轮按角度填充个体 并收集剩余个体
	{
		//initial_indi_pop[n].yobj2angle(ObserverPoint);
		initial_indi_pop[n].obj2angle_index( IdealPoint, popsize);

		if (sectorpopFlag[initial_indi_pop[n].sectorialindex] == false)//if (sectorialgrid[initial_indi_pop[n].sectorialindex] == -1)
		{
			population[initial_indi_pop[n].sectorialindex] = initial_indi_pop[n];
			sectorpopFlag[initial_indi_pop[n].sectorialindex] = true;
		}else{
			bool bReplaced;
			rest_ind = population[initial_indi_pop[n].sectorialindex];//.indiv;
			bReplaced = Pareto_HyperVolume_compare_sectorialgrid(initial_indi_pop[n]);
			if(bReplaced)
			{
				initial_rest_indi_pop.push_back( rest_ind );
			}else{
				initial_rest_indi_pop.push_back( initial_indi_pop[n] );
			}
		}
	}
	initial_indi_pop.clear();

	for(int n=0;n< popsize;n++)//填充剩余的个体
	{
		if (sectorpopFlag[ n ] == false)
		{
			double min_angle_difference, angle_difference ;
			int min_angle_difference_ind ;
			int rest_size = initial_rest_indi_pop.size();

			min_angle_difference = fabs( initial_rest_indi_pop[0].sectorialangle - (double)n / (double)(popsize-1) );
			//min_angle_difference = abs( initial_rest_indi_pop[0].sectorialindex - n );
			min_angle_difference_ind = 0 ;
			for(int j = 1; j<rest_size; j++ )
			{
				angle_difference = fabs( initial_rest_indi_pop[j].sectorialangle - (double)n / (double)(popsize-1) );
				//angle_difference = abs( initial_rest_indi_pop[j].sectorialindex - n );
				if(angle_difference < min_angle_difference)
				{
					min_angle_difference = angle_difference;
					min_angle_difference_ind = j;
				}
			}
			population[n] = initial_rest_indi_pop[min_angle_difference_ind];
			sectorpopFlag[n] = true;
			rest_ind = initial_rest_indi_pop[min_angle_difference_ind];
			initial_rest_indi_pop[min_angle_difference_ind] = initial_rest_indi_pop[ rest_size-1 ];
			initial_rest_indi_pop[ rest_size-1 ] = rest_ind;
			initial_rest_indi_pop.pop_back();
		}
	}

	//double min_angle_difference, angle_difference ;
	//int min_angle_difference_ind ;
	//int min_angle_insert_point;

	//for(int rest_size = initial_rest_indi_pop.size();rest_size>0;rest_size--)
	//{
	//	min_angle_difference = 2.0;
	//	for(int n=0;n< popsize;n++)//填充剩余的个体
	//	{
	//		if (sectorpopFlag[ n ] == false)
	//		{
	//			//min_angle_difference = abs( initial_rest_indi_pop[0].sectorialangle - (double)n / (double)(popsize-1) );
	//			//min_angle_difference_ind = 0 ;
	//			for(int j = 0; j<rest_size; j++ )
	//			{
	//				angle_difference = abs( initial_rest_indi_pop[j].sectorialangle - (double)n / (double)(popsize-1) );
	//				if(angle_difference < min_angle_difference)
	//				{
	//					min_angle_difference = angle_difference;
	//					min_angle_difference_ind = j;
	//					min_angle_insert_point = n;
	//				}
	//			}
	//		}
	//	}

	//	sectorpop[min_angle_insert_point] = initial_rest_indi_pop[min_angle_difference_ind];
	//	sectorpopFlag[min_angle_insert_point] = true;
	//	initial_rest_indi_pop[min_angle_difference_ind] = initial_rest_indi_pop[ rest_size-1 ];
	//}

	initial_rest_indi_pop.clear();
	sectorpopFlag.clear();
}

bool CCAEA::Pareto_HyperVolume_compare_sectorialgrid(CCAEAInd& ind)
{
	bool bReplaced = false;

	double contributing1, contributing2;

	//int ind_index = ind.sectorialindex;
	//TCompare iResult1 = ind.Compare( sectorpop[ind.sectorialindex] );//.indiv
	//if( sectorpop[ind.sectorialindex].angle_in_sector(ObserverPoint, popsize, anglesingle, ind.sectorialindex) )
	population[ind.sectorialindex].obj2angle_index(IdealPoint,popsize);
	if(population[ind.sectorialindex].sectorialindex == ind.sectorialindex )
	{
		//if(iResult1 == _Pareto_Dominating)
		//{
		//	sectorpop[ind.sectorialindex] = ind;//replace   //.indiv
		//	bReplaced = true;
		//}
		//else if(iResult1 == _Pareto_Nondominated)
		//{
			contributing1 = GetFastigiateHyperVolume(ind, ind.sectorialindex, ReferencePoint);
			contributing2 = GetFastigiateHyperVolume(population[ind.sectorialindex], ind.sectorialindex, ReferencePoint);

			if(contributing1 < contributing2)
			{
				population[ind.sectorialindex] = ind;//replace
				bReplaced = true;
			}
		//}
	}else{
		//if( iResult1 != _Pareto_Dominated )
		{
			population[ind.sectorialindex] = ind;//replace
			bReplaced = true;
		}
	}

	return bReplaced;
}

double CCAEA::GetFastigiateHyperVolume(CCAEAInd&  ind, int ind_index , vector <double> &ReferencePoint)
{
	double FastigiateVolume;

	double normalizedf[2];
	for(int j=0; j<nobj; j++)
	{
		normalizedf[j] = ind.y_obj[j] - IdealPoint[j];
	}
	if(ind_index==0){
		FastigiateVolume = 0.5*(ind_index+0.5)/(popsize-ind_index-1.5)*normalizedf[1]*normalizedf[1] + (ReferencePoint[1]-normalizedf[1])*normalizedf[0];
	}else if(ind_index == popsize-1){
		FastigiateVolume = 0.5*(popsize-ind_index-0.5)/(ind_index-0.5)*normalizedf[0]*normalizedf[0] + (ReferencePoint[0]-normalizedf[0])*normalizedf[1];
	}else{
		FastigiateVolume = 0.5*(popsize-ind_index-0.5)/(ind_index-0.5)*normalizedf[0]*normalizedf[0] + 0.5*(ind_index+0.5)/(popsize-ind_index-1.5)*normalizedf[1]*normalizedf[1] - normalizedf[0]*normalizedf[1];
	}

	return FastigiateVolume;
}

//double CCAEA::GetFastigiateHyperVolume(CCAEAInd&  ind, int ind_index, vector <double> &ReferencePoint)
//{
//	double FastigiateVolume;
//
//	double Volume1, Volume2, height1, height2, BottomArea1, BottomArea2;
//
//	double f1sharp, f2sharp;
//
//	double directionf1, directionf2;
//
//	height1 = ind.y_obj[1] - ObserverPoint[1];
//	if(ind_index == popsize-1 )
//	{
//		f1sharp = ReferencePoint[0];
//	}else{
//		directionf1 = (ind_index + 0.5) / (double) (popsize-1) ;
//		f1sharp = directionf1 / (1 - directionf1) * (ind.y_obj[1] - ObserverPoint[1]) +  ObserverPoint[0];
//	}
//	BottomArea1 = f1sharp - ind.y_obj[0];
//	Volume1 = 0.5 * BottomArea1 * height1;
//
//	height2 = ind.y_obj[0] - ObserverPoint[0];
//	if(ind_index == 0 )
//	{
//		f2sharp = ReferencePoint[1];
//	}else{
//		directionf2 = (ind_index - 0.5) / (double) (popsize-1) ;
//		f2sharp =  (1 - directionf2) / directionf2 * (ind.y_obj[0] - ObserverPoint[0]) +  ObserverPoint[1];
//	}
//	BottomArea2 = f2sharp - ind.y_obj[1];
//	Volume2 = 0.5 * BottomArea2 * height2;
//
//	FastigiateVolume = Volume1 + Volume2;
//
//	return FastigiateVolume;
//}


double  CCAEA::tour_selection_hv_difference(int p, vector <CCAEAInd>  &mypopulation)
{
	//int p1 = int(rnd_uni(&rnd_uni_init)*popsize);
	//int p2 = int(rnd_uni(&rnd_uni_init)*popsize);
	int num = 0;
	double hv1 = GetFastigiateHyperVolume(mypopulation[p], mypopulation[p].sectorialindex,ReferencePoint);
	double hv3,hv4,hv_difference=0;
	int neighbor1,neighbor2;
	if(p-1>=0)
	{
		neighbor1 = p-1;
	}else neighbor1 = p+2;
	hv3 = GetFastigiateHyperVolume(mypopulation[neighbor1], mypopulation[neighbor1].sectorialindex,ReferencePoint);
	hv_difference += hv3 ;
	num++;

	if(p+1<=popsize-1)
	{
		neighbor2 = p+1;
	}else neighbor2 = p-2;
	hv4 = GetFastigiateHyperVolume(mypopulation[neighbor2], mypopulation[neighbor2].sectorialindex,ReferencePoint);
	hv_difference += hv4 ;
	num++;

	hv_difference = hv_difference / num - hv1;

	return hv_difference;
}

int  CCAEA::tour_selection_hv2(vector <CCAEAInd>  &mypopulation)
{
	int p1 = int(rnd_uni(&rnd_uni_init)*(pops_sub_self[mpi_rank] -1) + pops_sub_start[mpi_rank]);
	int p2 = int(rnd_uni(&rnd_uni_init)*(pops_sub_self[mpi_rank] - (p1 - pops_sub_start[mpi_rank] +1) ) + p1 + 1);
	//if(p2>=p1) p2++;

	double hv1 = tour_selection_hv_difference( p1, mypopulation );
	double hv3 = tour_selection_hv_difference( p2, mypopulation );
	if(hv1 >= hv3 )
		return p1;
	else
		return p2;
	//double hv1 = GetFastigiateHyperVolume( mypopulation[p1], mypopulation[p1].sectorialindex, ReferencePoint );
	//double hv3 = GetFastigiateHyperVolume( mypopulation[p2], mypopulation[p2].sectorialindex, ReferencePoint );
	//if(hv1 <= hv3 )
	//	return p1;
	//else

	//	return p2;
}

void CCAEA::initTopologies()
{
	if (numprocs > 1)
	{
		//生成拓扑结构
		//单向环
		if (paraTopology == UNI_RING)
		{
			neighbor_partitions.push_back((mpi_rank + 1 + numprocs) % numprocs);
		}
		else if (paraTopology == BI_RING )
		{
			//双向环
			neighbor_partitions.push_back((mpi_rank + 1 + numprocs) % numprocs);
			neighbor_partitions.push_back((mpi_rank - 1 + numprocs) % numprocs);
		}
		else if (paraTopology == NEI)
		{	//
			if (mpi_rank == 0)
				neighbor_partitions.push_back((mpi_rank + 1 + numprocs) % numprocs);
			else if (mpi_rank == numprocs - 1)
				neighbor_partitions.push_back((mpi_rank - 1 + numprocs) % numprocs);
			else
			{
				neighbor_partitions.push_back((mpi_rank + 1 + numprocs) % numprocs);
				neighbor_partitions.push_back((mpi_rank - 1 + numprocs) % numprocs);
			}
		}
		else if (paraTopology == FULL)
		{
			//完全拓扑
			for (int i = 0; i < numprocs; i++)
			{
				if (mpi_rank != i)
					neighbor_partitions.push_back(i);
			}
		}
		partition_niche = neighbor_partitions.size();
	}
	/*
	cout << "["<<mpi_rank<<"] nei size : " << partition_niche << " : ";
	for (int i = 0; i < partition_niche; i++)
		cout << neighbor_partitions[i] << "  ";
	cout << endl;
	*/
}

void CCAEA::execute(int mg, int irun, vector<double>& igd, vector<double>& hv, double &runtime, PTYPE ptype)
{
	////seed = (seed + 23)%1377;
	////rnd_uni_init = -(long)seed;
	//vector<double> gd;
	initTopologies();
	if (numprocs >1) {
		//分配缓冲空间
		//更新理想点使用的缓冲区
		if (ptype == IDEAL_INDIV || ptype == IDEAL)
		{
			recv_idealpoint = new double*[partition_niche];
			send_idealpoint = new double*[partition_niche];
			req_idealpoint_recv = new MPI_Request[partition_niche];
			req_idealpoint_send = new MPI_Request[partition_niche];
			isBegin_Ideal = new bool[partition_niche];
		}
		//个体更新接收
		//决策向量+目标向量+权重向量
		s_indiv = 1 + nvar + nobj ;
		indivBufSize = 1 + s_indiv * popsize;//sizeof(int)+(sizeof(int) + sizeof(double)*(nvar+nobj)) * maxNiche;
		if (ptype == IDEAL_INDIV || ptype == INDIV)
		{
			recv_indiv = new double*[partition_niche];
			send_indiv = new double*[partition_niche];
			req_indiv_recv = new MPI_Request[partition_niche];
			req_indiv_send = new MPI_Request[partition_niche];
			isBegin_Indiv = new bool[partition_niche];
			indiv_buf = new double[indivBufSize];
		}
		//cout << "[" << mpi_rank << "]  3" << endl;
		//cout << "[" << mpi_rank << "] partition niche :" << partition_niche << endl;
		for (int i = 0; i < partition_niche; i++)
		{
			if (ptype == IDEAL_INDIV || ptype == IDEAL)
			{
				recv_idealpoint[i] = new double[nobj];
				send_idealpoint[i] = new double[nobj];
				isBegin_Ideal[i] = false;
				//非阻塞接收理想点更新
				MPI_Irecv(recv_idealpoint[i], nobj, MPI_DOUBLE, neighbor_partitions[i], TAG_IDEALPOINT, MPI_COMM_WORLD, &(req_idealpoint_recv[i]));
			}

			if (ptype == IDEAL_INDIV || ptype == INDIV)
			{
				recv_indiv[i] = new double[indivBufSize];
				send_indiv[i] = new double[indivBufSize];
				isBegin_Indiv[i] = false;
				MPI_Irecv(recv_indiv[i], indivBufSize, MPI_DOUBLE, neighbor_partitions[i], TAG_INDIV, MPI_COMM_WORLD, &(req_indiv_recv[i]));
			}
		}
		//cout<<"["<<mpi_rank<<"] "<<endl;
		MPI_Barrier(MPI_COMM_WORLD);
	}

	isIdealpointUpdate = false;
	recordUpdateofIndiv = new bool[population.size()];
	for (int i = 0; i < population.size(); i++)
		recordUpdateofIndiv[i] = false;

	init_population();

	if (ptype == IDEAL_INDIV || ptype == IDEAL)
	{
		if (numprocs > 1) {
			//cout << "[" <<mpi_rank << "] sync idealpoint" <<endl;
			for (int i = 0; i < nobj; i++)
			{
				send_idealpoint[0][i] = IdealPoint[i];
			}

			//cout << "[" << mpi_rank << "] ideal (" << idealpoint[0] << "," << idealpoint[1] << "," << idealpoint[2] << ")" << endl;

			double *recv_utopian;
			if (mpi_rank == 0)
				recv_utopian = new double[nobj*numprocs];

			//rank = 0 收集理想点
			//收集数据
			int *rcounts = new int[numprocs];
			for (int i = 0; i < numprocs; i++) {
				rcounts[i] = nobj;
			}
			int *displs = new int[numprocs];
			displs[0] = 0;
			for (int i = 1; i < numprocs; i++)
				displs[i] = displs[i - 1] + rcounts[i - 1];

			MPI_Gatherv(send_idealpoint[0], nobj, MPI_DOUBLE, recv_utopian, rcounts, displs, MPI_DOUBLE, root_rank, MPI_COMM_WORLD);

			//统一理想点
			if (mpi_rank == root_rank)
			{
				for (int i = 0; i < numprocs; i++)
				{
					int j;
					for (j = 0; j < nobj; j++)
					{
						double tmp = recv_utopian[i*nobj + j];
						if (tmp < IdealPoint[j])
							IdealPoint[j] = tmp;
					}
				}
				//cout << "tongyi(" << idealpoint[0] << "," << idealpoint[1] << "," << idealpoint[2] << ")" << endl;
			}
			//分发、同步理想点
			for (int i = 0; i < nobj; i++)
			{
				send_idealpoint[0][i] = IdealPoint[i];
			}
			MPI_Bcast(send_idealpoint[0], nobj, MPI_DOUBLE, root_rank, MPI_COMM_WORLD);

			//更新进程内理想点
			for (int i = 0; i < nobj; i++) {
				IdealPoint[i] = send_idealpoint[0][i];
			}
			//cout << "[" << mpi_rank << "] after sync idealpoint("<<idealpoint[0]<<","<<idealpoint[1]<<","<<idealpoint[2]<<")" << endl;
			if (mpi_rank == root_rank)
				delete[]recv_utopian;
			delete[]rcounts;
			delete[]displs;
		}
	}


	int gen = 1;
	double igd_gen = 0;
	double hv_gen = 0;
	char filename[1024];
	clock_t begintime;
	clock_t endtime;
	unevolvetime = 0;

	if (IS_DEBUG) {
		begintime = clock();
		//收集计算解集质量的必要数据,目标向量
		//等待，同步
		if (mpi_rank == 0)
		{
			//载入真实前沿
			sprintf(filename, "%s/TestData/PF_Real/%s(%d).dat", WORKINGDIR, strTestInstance, nobj);
			loadpfront(filename, ps);
		}

		if (MODE_DEBUG)
		{
			//cout << "[" << mpi_rank << "] before gather pop y" << endl;
			gather_pop_y();
			//cout << "[" << mpi_rank << "] after" << endl;
			//Root进程处理数据
			if (mpi_rank == root_rank)
			{
				//计算解集质量（IGD，HV)
				calc_qualities(igd_gen, hv_gen);
				igd.push_back(gen);  igd.push_back(igd_gen);
				hv.push_back(gen);  hv.push_back(hv_gen);
				cout << "gen =  " << gen << "  hypervolume = " << hv_gen << "  " << "  IGD = " << igd_gen << "  " << endl;
			}
			//MPI_Barrier(MPI_COMM_WORLD);
		}
		//等待，同步
		MPI_Barrier(MPI_COMM_WORLD);
		endtime = clock();
		unevolvetime += endtime - begintime;
	}

	int iLen = 0;
	int Len = ceil(popsize / 10.0);
	CCAEAInd child2;
	//开始执行进化
	for (gen = 2; gen <= mg; gen++)
	{
		for (int i = 0; i < pops_sub_self[mpi_rank]; i++)
		{
			int parent_index1, parent_index2;

			iLen = iLen % Len;
			if (iLen == 0)
				parent_index1 = 0;
			else if (iLen == 1)
				parent_index1 = popsize - 1;
			else
				parent_index1 = tour_selection_hv2(population);// int(rnd_uni(&rnd_uni_init)*( popsize )); ////int(rnd_uni(&rnd_uni_init)*( popsize )); //
			//ind_selected = &(population[parent_index1]);

			parent_index2 = tour_selection_hv2(population);// int(rnd_uni(&rnd_uni_init)*( popsize ));// //int(rnd_uni(&rnd_uni_init)*( popsize )); //
			//ind_selected2 = &(population[parent_index2]);

			//int parent_index3 =  tour_selection_hv2( population );

		//sbx_mutation_evolution(*ind_selected, *ind_selected2);
			//if(rnd_uni(&rnd_uni_init)<0.9)
			real_sbx_xoverA(population[parent_index1], population[parent_index2], onechild, child2);
			//diff_evo_xoverB(population[parent_index1], population[parent_index2], population[parent_index3], onechild, 0.5, 1.0);
			realmutation(onechild, 1.0 / nvar);

			onechild.obj_eval();

			//nfes++;
			iLen++;

			bool updated = update_extreme_point(onechild);//update_anglerange_maxminindi(onechild, true); //update_nondominated_obj_maxmin(onechild);//update_angle_maxmin(onechild);//
			if (updated)
			{
				reset_angle();
			}

			//onechild.yobj2angle(ObserverPoint);
			onechild.obj2angle_index(IdealPoint, popsize);//之前之所以出现奇异点 是因为计算onechild角度是在统计极端点前 极端点改变了 而onechild的角度是老的

			//update_sectorialgrid(onechild);
			if (Pareto_HyperVolume_compare_sectorialgrid(onechild))
				recordUpdateofIndiv[onechild.sectorialindex] = true;
		}



		if (ptype == IDEAL_INDIV || ptype == IDEAL)
		{
			if (numprocs > 1)
			{
				cout <<gen<< "    [" << mpi_rank << "] update idealpoint" << endl;
				if (isIdealpointUpdate)
				{
					//cout << "has update" << endl;
					for (int i = 0; i < partition_niche; i++)
					{
						sendIdealpoint(send_idealpoint[i], isBegin_Ideal[i], neighbor_partitions[i], &(req_idealpoint_send[i]));
					}
					isIdealpointUpdate = false;
				}
				//cout << "[" << mpi_rank << "] update idealpoint finish" << endl;
			}
			update_idealpoint_from_outside();
		}

		if (ptype == IDEAL_INDIV || ptype == INDIV)
		{
			if (numprocs > 1 && gen % comm_i == 0)
			{
				cout << "[" << mpi_rank << "] update indiv" << endl;
				//更新个体

				int updated_size = 0;
				int m = 1;
				for (int i = 0; i < population.size(); i++)
				{
					if ( (i < pops_sub_start[mpi_rank] || i > (pops_sub_start[mpi_rank] + pops_sub_self[mpi_rank]) ) && (recordUpdateofIndiv[i] == true))
					{
						indiv_buf[m++] = i;
						for (int v = 0; v < nvar; v++)
							indiv_buf[m++] = population[i].x_var[v];
						for (int v = 0; v < nobj; v++)
							indiv_buf[m++] = population[i].y_obj[v];
						updated_size++;
					}
				}

				if (updated_size > 0)
				{
					indiv_buf[0] = updated_size;
					int send_size = 1 + updated_size * s_indiv;
					for (int i = 0; i < partition_niche; i++)
					{
						int neiProc = neighbor_partitions[i];

						//memcpy(send_indiv[i],indiv_buf, send_size * sizeof(double));
						for (int k = 0; k < 1 + updated_size*(s_indiv); k++)
							send_indiv[i][k] = indiv_buf[k];
						/*
						cout << "[" << mpi_rank << "] up size:" << send_indiv[i][0] << "     ";
						for (int k = 1; k < send_indiv[i][0]*(1 + nvar + nobj); k++)
						cout << send_indiv[i][k] << ",";
						cout << endl;
						*/
						if (isBegin_Indiv[i])
						{
							//cout << gen <<"   [" << mpi_rank << "] wait indiv" << endl;
							MPI_Wait(&req_indiv_send[i], &status_indiv);
							//MPI_Cancel(&req_indiv_send[i]);
							//cout << gen<<"   [" << mpi_rank << "] wait indiv finish" << endl;
						}
						//cout << gen<< "   [" << mpi_rank << "] -->" << neiProc << "   send indiv" << endl;
						MPI_Isend(send_indiv[i], send_size, MPI_DOUBLE, neiProc, TAG_INDIV, MPI_COMM_WORLD, &(req_indiv_send[i]));
						isBegin_Indiv[i] = true;
					}
				}

				//cout<<"["<<mpi_rank<<"] 5"<<endl;

				//重置更新标志
				for (int i = 0; i < population.size(); i++)
					recordUpdateofIndiv[i] = false;

			}
			update_indiv_from_outside();
		}

		if (IS_DEBUG && gen % I_OUT == 0) {

			//获取计算解集质量的必要数据，目标向量
			//等待，同步
			//update_problem_from_outside();
			//MPI_Barrier(MPI_COMM_WORLD);
			//update_problem_from_outside();
			if (MODE_DEBUG)
			{
				//等待，同步
				//MPI_Barrier(MPI_COMM_WORLD);
				begintime = clock();
				cout << gen<<"   [" << mpi_rank << "] before gatherv" << endl;
				gather_pop_y();
				cout << "[" << mpi_rank << "] after gather" << endl;
				//Root进程处理数据
				if (mpi_rank == root_rank)
				{
					//计算解集质量（IGD，HV)
					calc_qualities(igd_gen, hv_gen);
					igd.push_back(gen);  igd.push_back(igd_gen);
					hv.push_back(gen);  hv.push_back(hv_gen);
					cout << "gen =  " << gen << "  hypervolume = " << hv_gen << "  " << "  IGD = " << igd_gen << "  " << endl;
				}
				//等待，同步
				MPI_Barrier(MPI_COMM_WORLD);
				endtime = clock();
				unevolvetime += endtime - begintime;
			}

		}

	}

	if (numprocs > 1)
	{
		if (ptype == IDEAL_INDIV || ptype == INDIV)
		{
			isDoneRecv = new int[partition_niche];
			isDoneSend = new int[partition_niche];
			b_isdone = new int[partition_niche];
			req_done_recv = new MPI_Request[partition_niche];
			req_done_send = new MPI_Request[partition_niche];

			//清理操作，完成操作
			for (int i = 0; i < partition_niche; i++)
			{
				isDoneRecv[i] = 0;
				int flag;
				MPI_Irecv(&b_isdone[i], 1, MPI_INT, neighbor_partitions[i], TAG_DONE, MPI_COMM_WORLD, &req_done_recv[i]);
				if (isBegin_Indiv[i])
				{
					MPI_Test(&req_indiv_send[i], &flag, &status_indiv);
					if (flag)
						isBegin_Indiv[i] = false;
				}
				if (isBegin_Indiv[i])
					isDoneSend[i] = 0;
				else
				{
					isDoneSend[i] = 1;
					MPI_Isend(&isDoneSend[i], 1, MPI_INT, neighbor_partitions[i], TAG_DONE, MPI_COMM_WORLD, &req_done_send[i]);
				}
			}


			while (1)
			{
				for (int i = 0; i < partition_niche; i++)
				{
					int flag;
					if (!isDoneSend[i])
					{
						if (isBegin_Indiv[i])
						{
							MPI_Test(&req_indiv_send[i], &flag, &status_indiv);
							if (flag)
								isBegin_Indiv[i] = false;
						}
						if (!isBegin_Indiv[i])
						{
							isDoneSend[i] = 1;
							MPI_Isend(&isDoneSend[i], 1, MPI_INT, neighbor_partitions[i], TAG_DONE, MPI_COMM_WORLD, &req_done_send[i]);
						}
					}
					if (!isDoneRecv[i])
					{
						MPI_Test(&req_done_recv[i], &flag, &status_done);
						if (flag && b_isdone[i])
							isDoneRecv[i] = 1;
						else
						{
							if (ptype == IDEAL_INDIV)
								update_idealpoint_from_outside();

							getAndUpdateIndivs(recv_indiv[i], indivBufSize, neighbor_partitions[i], &req_indiv_recv[i]);
						}
					}

				}
				bool isFinish = true;
				for (int i = 0; i < partition_niche; i++)
				{
					if (!isDoneRecv[i])
					{
						isFinish = false;
						break;
					}
					if (!isDoneSend[i])
					{
						isFinish = false;
						break;
					}
				}
				if (isFinish)
					break;
			}
		}

		MPI_Barrier(MPI_COMM_WORLD);

		if (ptype == IDEAL_INDIV || ptype == IDEAL)
			update_idealpoint_from_outside();

		for (int i = 0; i < partition_niche; i++)
		{
			if (ptype == IDEAL_INDIV || ptype == IDEAL)
			{
				MPI_Cancel(&req_idealpoint_recv[i]);
				MPI_Request_free(&req_idealpoint_recv[i]);
				//等待，同步
				if (isBegin_Ideal[i])
				{
					int flag = 1;
					MPI_Test(&req_idealpoint_send[i], &flag, &status_idealpoint);
					if (flag == 0)
					{
						MPI_Cancel(&req_idealpoint_send[i]);
						MPI_Request_free(&req_idealpoint_send[i]);
					}
				}
			}

			if (ptype == IDEAL_INDIV || ptype == INDIV)
			{
				MPI_Cancel(&req_indiv_recv[i]);
				MPI_Request_free(&req_indiv_recv[i]);
				if (isBegin_Indiv[i])
				{
					int flag;
					MPI_Test(&req_indiv_send[i], &flag, &status_indiv);
					if (flag)
					{
						MPI_Cancel(&req_indiv_send[i]);
						MPI_Request_free(&req_indiv_send[i]);
					}
				}
			}
		}

	}

	MPI_Barrier(MPI_COMM_WORLD);
	finish = clock();
	runtime = (double)(finish - start - unevolvetime) / CLOCKS_PER_SEC;

	//cout << "[" << mpi_rank << "] idealpoint("<<IdealPoint[0]<<","<<IdealPoint[1]<<")" << endl;
	//cout << "[" << mpi_rank << "] TrueNadirpoint(" << TrueNadirPoint[0] << "," << TrueNadirPoint[1] << ")" << endl;
	//cout << "[" << mpi_rank << "] Referentpoint(" << ReferencePoint[0] << "," << ReferencePoint[1] << ")" << endl;


	if (IS_DEBUG) {
		if (!MODE_DEBUG || (MODE_DEBUG && (gen - 1) % I_OUT != 0))
		{//最后一代的数据
		 //获取计算解集质量的必要数据，目标向量
		 //等待，同步
		 //MPI_Barrier(MPI_COMM_WORLD);
			gather_pop_y();
			//Root进程处理数据
			if (mpi_rank == 0)
			{
				//计算解集质量（IGD，HV)
				calc_qualities(igd_gen, hv_gen);
				igd.push_back(gen - 1);  igd.push_back(igd_gen);
				hv.push_back(gen - 1);  hv.push_back(hv_gen);
				//cout<<"gen =  "<<gen-1<<"  hypervolume = "<<hv_gen<<"  "<<"  IGD = "<<igd_gen<<"  "<<endl;
			}
		}
		/*
		// save the final population - X space
		//每个进程都保存，用于后面展示
		sprintf(filename, "%s/TestData/%d/POP/POP_%s_%s(%d)_%d_%dR%d.dat",WORKINGDIR,mpi_rank,STR_MOEAD,strTestInstance,nobj, popsize, total_gen, irun);
		save_front(filename);
		sprintf(filename, "%s/TestData/%d/POF/PF_%s_%s(%d)_%d_%dR%d.dat",WORKINGDIR,mpi_rank,STR_MOEAD, strTestInstance,nobj, popsize, total_gen, irun);
		vector <CIndividual>  frontpop;
		population2front(population, frontpop);
		save_population(frontpop, filename);//保存扇形存档种群的非劣前沿

		frontpop.clear();

		sprintf(filename, "%s/TestData/%d/POS/POS_%s_%s(%d)_%d_%dR%d.dat",WORKINGDIR,mpi_rank,STR_MOEAD, strTestInstance,nobj, popsize, total_gen, irun);
		save_pos(filename);
		*/
	}
	//输出结果
	//合并前沿数据
	//等待，同步
	//	MPI_Barrier(MPI_COMM_WORLD);
	//cout << "[" << mpi_rank << "] gather" << endl;
	gather_populations();
	//cout << "[" << mpi_rank << "] gather end" << endl;
	//等待，同步
	//cout<<"["<<mpi_rank<<"] send data"<<endl;
	//MPI_Barrier(MPI_COMM_WORLD);
	if (mpi_rank == 0)
	{
		//cout << "receive mark" << endl;
		//cout<<"["<<mpi_rank<<"] 处理数据"<<endl;
		//处理数据
		vector<vector<double> > pop_x;
		vector<vector<double> > pop_y;
		vector<int> pop_index;
		pop_x.reserve(popsize);
		pop_y.reserve(popsize);
		int m = 0;
		for (int i = 0; i<popsize; i++)
		{
			vector<double> x_var;
			for (int j = 0; j<nvar; j++)
				x_var.push_back(b_recv[m++]);
			pop_x.push_back(x_var);
			vector<double> y_obj;
			for (int j = 0; j<nobj; j++)
				y_obj.push_back(b_recv[m++]);
			pop_y.push_back(y_obj);
			pop_index.push_back(i);
			//cout << "(" << y_obj[0] << "," << y_obj[1] << ")" << endl;
		}

		//cout << "receive end" << endl;
		vector<int> index_split;
		for (int i = 0; i<numprocs; i++)
		{
			for (int j = 0; j<pops_sub_self[i]; j++)
				index_split.push_back(i);
		}

		//cout << "dewfrwev" << endl;
		//保存种群前沿，未清理
		sprintf(filename, "%s/TestData/POP/POP_CAEA_%s(%d)_%d_%dR%d.dat", WORKINGDIR, strTestInstance, nobj, popsize, total_gen, irun);
		save_population_front_split(pop_y, pop_index, index_split, filename);//保存扇形存档种群的非劣前沿
		pop_index.clear();
		//过滤取出前沿
		vector <vector<double> >  pf_y;
		vector<int> pf_index;
		population2front(pop_y, pf_y, pf_index);
		//保存整体前沿
		sprintf(filename, "%s/TestData/POF/PF_CAEA_%s(%d)_%d_%dR%d.dat", WORKINGDIR, strTestInstance, nobj, popsize, total_gen, irun);
		save_population_front_split(pf_y, pf_index, index_split, filename);//保存扇形存档种群的非劣前沿
																		   //保存决策变量
		sprintf(filename, "%s/TestData/POS/POS_CAEA_%s(%d)_%d_%dR%d.dat", WORKINGDIR, strTestInstance, nobj, popsize, total_gen, irun);
		save_pos_split(pop_x, pf_index, index_split, filename);//保存扇形存档种群的非劣前沿

		pop_y.clear();
		pf_y.clear();
		index_split.clear();
		pf_index.clear();
	}

	population.clear();

	if (numprocs > 1)
	{
		//cout <<"["<<mpi_rank<<"] finall" << endl;
		for (int i = 0; i < partition_niche; i++)
		{
			if (ptype == IDEAL_INDIV || ptype == IDEAL)
			{
				delete[]recv_idealpoint[i];
				recv_idealpoint[i] = NULL;
				delete[]send_idealpoint[i];
				send_idealpoint[i] = NULL;
			}
			if (ptype == IDEAL_INDIV || ptype == INDIV)
			{
				delete[]recv_indiv[i];
				recv_indiv[i] = NULL;
				delete[]send_indiv[i];
				send_indiv[i] = NULL;
				delete[]indiv_buf;
				indiv_buf = NULL;
			}
		}
		if (ptype == IDEAL_INDIV || ptype == IDEAL)
		{
			delete[]recv_idealpoint;
			recv_idealpoint = NULL;
			delete[]send_idealpoint;
			send_idealpoint = NULL;
			delete[]req_idealpoint_recv;
			req_idealpoint_recv = NULL;
			delete[]req_idealpoint_send;
			req_idealpoint_send = NULL;
			delete[]isBegin_Ideal;
			isBegin_Ideal = NULL;
		}
		if (ptype == IDEAL_INDIV || ptype == INDIV)
		{
			delete[]req_indiv_recv;
			req_indiv_recv = NULL;
			delete[]req_indiv_send;
			req_indiv_send = NULL;
			delete[]isBegin_Indiv;
			isBegin_Indiv = NULL;

			delete[]req_done_recv;
			req_done_recv = NULL;
			delete[]req_done_send;
			req_done_send = NULL;
			delete[]isDoneRecv;
			isDoneRecv = NULL;
			delete[]isDoneSend;
			isDoneSend = NULL;
			delete[]b_isdone;
			b_isdone = NULL;
		}

		neighbor_partitions.clear();
		//cout << "[" << mpi_rank << "]finish" << endl;
	}
}

//void CCAEA::uniform_selection(CCAEAInd*& ind_selected)
//{
//	int pop_index;
//
//	pop_index = int(rnd_uni(&rnd_uni_init)*( popsize ));
//
//	ind_selected = &(population[pop_index]);//.indiv
//}


void CCAEA::reset_angle()
{
	int n;
	for(n=0; n<popsize; n++)
	{
		//population[n].yobj2angle(ObserverPoint);
		population[n].obj2angle_index(IdealPoint, popsize);
	}
}


bool CCAEA::update_idealpoint_from_outside()
{
	bool isUpdate = false;
	if (numprocs >1) {
		//只有多进程才使用，保证单进程可以执行
		for (int i = 0; i < partition_niche; i++)
		{
			int flag = 0;
			bool tmpUpdate = false;
			MPI_Test(&(req_idealpoint_recv[i]), &flag, &status_idealpoint);
			while (flag != 0)
			{//接收到来自其他进程的理想点更新

			 //cout << "[" << mpi_rank << "]("<<recv_idealpoint[i][0]<<","<<recv_idealpoint[i][1]<<","<<recv_idealpoint[i][2]<<")" << endl;
				tmpUpdate = update_ideal_point(recv_idealpoint[i]);
				//cout<<"["<<mpi_rank<<"]("<<recv_idealpoint[0]<<","<<recv_idealpoint[1]<<") old("<<ij[0]<<","<<ij[1]<<")now("<<idealpoint[0]<<","<<idealpoint[1]<<")"<<endl;

				//继续等待接收更新
				MPI_Irecv(recv_idealpoint[i], nobj, MPI_DOUBLE, neighbor_partitions[i], TAG_IDEALPOINT, MPI_COMM_WORLD, &(req_idealpoint_recv[i]));
				flag = 0;
				MPI_Test(&(req_idealpoint_recv[i]), &flag, &status_idealpoint);
			}
			isUpdate |= tmpUpdate;
		}
	}
	return isUpdate;
}

void CCAEA::sendIdealpoint(double *buf, bool &isSend, int target, MPI_Request* request)
{
	if (numprocs >1)
	{
		//多进程才使用，保证单进程也可以执行
		//如果理想点更新，则向下一个进程发送更新信息
		//阻塞，等待上一个理想点更新信息传送完成

		//可以发送更新信息
		for (int i = 0; i<nobj; i++)
		{
			buf[i] = IdealPoint[i];
		}
		MPI_Isend(buf, nobj, MPI_DOUBLE, target, TAG_IDEALPOINT, MPI_COMM_WORLD, request);
		//发送指令完毕，非阻塞，继续执行接下来计算,边传输边计算
		//cout<<"["<<mpi_rank<<"] send"<<endl;
		isSend = true;
	}
}


void CCAEA::getAndUpdateIndivs(double* b_recv_indiv, int datasize, int source, MPI_Request *request)
{
	//if(mpi_rank == 0)
	//cout<<"["<<mpi_rank<<"] LSize : "<<LSize<<"\tRSize : "<<RSize<<endl;
	int flag = 0;
	MPI_Test(request, &flag, &status_indiv);
	if (flag != 0)
		//if(flag !=0)
	{//接收到来自进程的个体更新
	 //从缓存中获取数据
		int m = 0;
		int count = (int)b_recv_indiv[m++];
		for (int p = 0; p < count; p++)
		{
			CCAEAInd indiv;
			int index = (int)(b_recv_indiv[m++] + 0.5);
			for (int i = 0; i < nvar; i++)
				indiv.x_var[i] = b_recv_indiv[m++];
			for (int i = 0; i < nobj; i++)
				indiv.y_obj[i] = b_recv_indiv[m++];

			//update extreme point
			update_extreme_point(indiv);

			indiv.obj2angle_index(IdealPoint,popsize);

			//根据权重向量计算对应种群下标
			//cout<<"["<<mpi_rank<<"]"<<index<<endl;
			if (Pareto_HyperVolume_compare_sectorialgrid(indiv))
				recordUpdateofIndiv[indiv.sectorialindex] = false;
		}

		MPI_Irecv(b_recv_indiv, datasize, MPI_DOUBLE, source, TAG_INDIV, MPI_COMM_WORLD, request);
	}
}

void CCAEA::update_indiv_from_outside()
{
	if (numprocs > 1)
	{
		for (int i = 0; i < partition_niche; i++)
			getAndUpdateIndivs(recv_indiv[i], indivBufSize, neighbor_partitions[i], &(req_indiv_recv[i]));
	}

}


void CCAEA::save_front(char saveFilename[1024])
{
	std::fstream fout;
	fout.open(saveFilename, std::ios::out);
	for (int n = 0; n<population.size(); n++)
	{
		for (int k = 0; k<nobj; k++)
			fout << population[n].y_obj[k] << "  ";
		fout << "\n";
	}
	fout.close();
}

void CCAEA::save_pos(char saveFilename[1024])
{
	std::fstream fout;
	fout.open(saveFilename, std::ios::out);
	for (int n = 0; n<population.size(); n++)
	{
		for (int k = 0; k<nvar; k++)
			fout << population[n].x_var[k] << "  ";
		fout << "\n";
	}
	fout.close();
}

void CCAEA::operator=(const CCAEA &emo)
{
	popsize = emo.popsize;
	population = emo.population;
}

void CCAEA::gather_pop_y()
{
	//发送各个个体的目标向量，用于计算IGD和HV指标
	int m = 0;
	int endPos = pops_sub_start[mpi_rank] + pops_sub_self[mpi_rank] ;
	for (int i = pops_sub_start[mpi_rank] ; i< endPos; i++)
	{
		//for(int j =0;j<nvar;j++)
		//b_send[m++] = population[i].indiv.x_var[j];
		for (int j = 0; j<nobj; j++) {
			b_send[m++] = population[i].y_obj[j];
		}
	}
	int size_send = nobj * pops_sub_self[mpi_rank];
	//cout << "[" << mpi_rank << "] m : " << m << "  count : " << size_send << endl;
	/*
	for (int i = 0; i < m; i++)
	{
		if (i % 2 == 0) cout << ",("<<i/2<<")";
		cout << b_send[i] << " ";
	} */
	//cout << endl;

	//收集数据
	/*int *rcounts = new int[numprocs];
	for (int i = 0; i<numprocs; i++)
		rcounts[i] = nobj*pops_sub_self[i];
	int *displs = new int[numprocs];
	displs[0] = 0;
	for (int i = 1; i<numprocs; i++)
		displs[i] = displs[i - 1] + rcounts[i - 1];
			 /*
	for (int i = 0; i < 2 * popsize; i++)
	{
		cout <<"("<<i<<")"<< b_recv[i] << " ";
	}	   */
	//cout << endl;
	/*
	cout <<root_rank<< "[" << mpi_rank << "] rcount : ";
	for (int i = 0; i < numprocs; i++)
		cout << rcounts[i] << " ";
	cout << endl;
	cout << "[" << mpi_rank << "] displs : ";
	for (int i = 0; i < numprocs; i++)
		cout << displs[i] << " ";
	cout << endl;
	  */
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Gatherv(b_send, size_send, MPI_DOUBLE, b_recv, rcounts_y, displs_y, MPI_DOUBLE, root_rank, MPI_COMM_WORLD);
	//cout << "after " << endl;
			   /*
	delete rcounts;
	rcounts = NULL;
	delete displs;
	displs = NULL;
				*/
	//cout << "out" << endl;
}

void CCAEA::gather_populations()
{
	//cout << "["<<mpi_rank<<"]mark" << endl;
	//发送各个个体的目标向量，用于计算IGD和HV指标
	int m = 0;
	int endPos = pops_sub_start[mpi_rank] + pops_sub_self[mpi_rank];
	for (int i = pops_sub_start[mpi_rank]; i < endPos; i++)
	{
		for (int j = 0; j < nvar; j++)
			b_send[m++] = population[i].x_var[j];
		for (int j = 0; j < nobj; j++) {
			b_send[m++] = population[i].y_obj[j];
		}
	}
	//cout <<"["<<mpi_rank<<"]"<<s_send<< "  adasd:"<<LLSize << "\t" << LLSize + pops_sub_self[mpi_rank] << endl;
	//收集数据
	/*
	int *rcounts = new int[numprocs];
	for (int i = 0; i<numprocs; i++)
		rcounts[i] = (nvar + nobj)*pops_sub_self[i];
	int *displs = new int[numprocs];
	displs[0] = 0;
	for (int i = 1; i<numprocs; i++)
		displs[i] = displs[i - 1] + rcounts[i - 1];
	*/
	MPI_Gatherv(b_send, s_send, MPI_DOUBLE, b_recv, rcounts_a, displs_a, MPI_DOUBLE, root_rank, MPI_COMM_WORLD);
	/*
	delete rcounts;
	rcounts = NULL;
	delete displs;
	displs = NULL;
	*/
	//cout << "["<<mpi_rank<<"]mark end" << endl;
}

void CCAEA::calc_qualities(double &igd, double& hv)
{
	//只有Root进程才使用
	//处理数据
	vector<vector<double> > pop_y;
	pop_y.reserve(popsize);
	int m = 0;
	for (int i = 0; i<popsize; i++)
	{
		vector<double> y_obj;
		cout << "dasda" <<endl;
		for (int j = 0; j<nobj; j++)
			y_obj.push_back(b_recv[m++]);
		cout << "[" << mpi_rank << "](" << y_obj[0] << "," << y_obj[1]  << ")" << endl;;
		pop_y.push_back(y_obj);
	}

	//过滤取出前沿
	vector <vector<double> >  pf_y;
	population2front(pop_y, pf_y);
	//计算IGD
	igd = calc_distance(pf_y);
	//计算HV
	hv = calc_hypervolume(pf_y);
	pop_y.clear();
	pf_y.clear();
	/*//只有Root进程才使用
 //处理数据
	cout << "[" << mpi_rank << "] begin cal" << endl;
	vector< vector<double> > pop_y;// (popsize);
	cout << "[" << mpi_rank << "] begin cal1" << endl;
	//pop_y.reserve(popsize);
	int m = 0;
	int vbn = 0;

	for (int i = 0; i < popsize; i++)
	{
		cout << "[" << mpi_rank << "] count : " << i <<"  obj : "<<nobj<< endl;
		vector<double> y_obj;
		for (int j = 0; j < nobj; j++)
		{
			double tmp = b_recv[m++];
			cout << "recv ("<<j<<"): " << tmp << endl;
			y_obj.push_back(tmp);

		}

		cout << "[" << mpi_rank << "](" << y_obj[0] << "," << y_obj[1] << ")" << endl;
		//pop_y[vbn] = y_obj;
		pop_y.push_back(y_obj);
		y_obj.clear();
		vbn++;
	}
	cout << "[" << mpi_rank << "] m:"<<m<<"vbn:" << vbn << endl;
	/*
	for (int i = 0; i < vbn; i++)
	{
		cout <<"("<<i<<")" <<pop_y[i][0] << " " << pop_y[i][1] << " , ";
	}

	cout << endl;

	//过滤取出前沿
	vector <vector<double> >  pf_y;
	cout << "[" << mpi_rank << "]1:"<<endl;
	population2front(pop_y, pf_y);

	cout << "[" << mpi_rank << "]2:" << endl;
	//计算IGD
	igd = calc_distance(pf_y);

	cout << "[" << mpi_rank << "]3:" << endl;
	//计算HV
	hv = calc_hypervolume(pf_y);

	cout << "[" << mpi_rank << "]4:" << endl;
	cout << "[" << mpi_rank << "] igd : "<<igd<<"    hv : "<<hv << endl;
	pop_y.clear();
	pf_y.clear();
*/
}

double CCAEA::calc_distance(vector<vector<double> >& ParetoFront)
{
	double distance = 0;
	for (int i = 0; i<ps.size(); i++)
	{
		double min_d = 1.0e+10;
		for (int j = 0; j<ParetoFront.size(); j++)
		{
			double d = dist_vector(ps[i].y_obj, ParetoFront[j]);
			if (d<min_d)  min_d = d;
		}
		distance += min_d;
	}
	distance /= ps.size();

	return distance;
}

double CCAEA::calc_hypervolume(vector<vector<double> >& ParetoFront)
{
	wfg hv;
	cout << "chechkjhlhjl" << endl;
	for (int i = 0; i < ParetoFront.size(); i++)
	{
		cout << ParetoFront[i][0] << "\t" << ParetoFront[i][1] << endl;
	}
	double result = hv.compute(ParetoFront, hv_ref_ub);

	cout << "dasdasd" << endl;
	return result / cubeHV;
}





void population2front(vector <CCAEAInd>  &mypopulation, vector <CCAEAInd>  &population_front)
{
	vector<int> nDominated;
	for (int n = 0; n<mypopulation.size(); n++)
		nDominated.push_back(0);


	for (int k = 0; k<mypopulation.size(); k++)
		for (int j = k + 1; j<mypopulation.size(); j++)
		{
			TCompare tresult = ParetoCompare(mypopulation[k].y_obj, mypopulation[j].y_obj);
			if (tresult == _Pareto_Dominated)
				nDominated[k]++;
			else if (tresult == _Pareto_Dominating)
				nDominated[j]++;
		}

	for (int n = 0; n<mypopulation.size(); n++)
		if (nDominated[n] == 0)
			population_front.push_back(mypopulation[n]);

	nDominated.clear();
}
void population2front(vector <CCAEAInd>  &mypopulation, vector <vector<double> >  &population_front)
{
	vector<int> nDominated;
	for (int n = 0; n<mypopulation.size(); n++)
		nDominated.push_back(0);


	for (int k = 0; k<mypopulation.size(); k++)
		for (int j = k + 1; j<mypopulation.size(); j++)
		{
			TCompare tresult = ParetoCompare(mypopulation[k].y_obj, mypopulation[j].y_obj);
			if (tresult == _Pareto_Dominated)
				nDominated[k]++;
			else if (tresult == _Pareto_Dominating)
				nDominated[j]++;
		}

	for (int n = 0; n<mypopulation.size(); n++)
		if (nDominated[n] == 0)
			population_front.push_back(mypopulation[n].y_obj);

	nDominated.clear();
}

void save_population(vector <CCAEAInd>  &mypopulation, char saveFilename[1024])
{
	std::fstream fout;
	fout.open(saveFilename, std::ios::out);
	for (int n = 0; n<mypopulation.size(); n++)
	{
		for (int k = 0; k<nobj; k++)
		{
			fout << mypopulation[n].y_obj[k] << "  ";
		}
		fout << "\n";
	}
	fout.close();
}

//使用图进行初始化数据
void CCAEA::init_population(MGraph *G)
{
	vector <bool>  sectorpopFlag;

	vector <CCAEAInd>  initial_indi_pop;
	CCAEAInd rest_ind;
	CCAEAInd ind;

	for (int n = 0; n<popsize; n++)//统计极值点or观测点  初始个体数目设置
	{
		sectorpopFlag.push_back(false);
		ind.rnd_init(G);
		//nfes++;
		initial_indi_pop.push_back(ind);//insert_sectorialgrid(ind, n);
	}

	//initial_observation_and_reference_point(initial_indi_pop[0]);

	for (int j = 0; j<nobj; j++)
	{
		IdealPoint.push_back(initial_indi_pop[0].y_obj[j]);
		TrueNadirPoint.push_back(initial_indi_pop[0].y_obj[j]);
		ReferencePoint.push_back(initial_indi_pop[0].y_obj[j]);
	}


	for (int n = 1; n<popsize; n++)//统计极值点or观测点  初始个体数目设置
	{
		update_extreme_point(initial_indi_pop[n]);
	}

	for (int n = 0; n<popsize; n++)//统计极值点or观测点  初始个体数目设置
	{
		population.push_back(initial_indi_pop[0]);
	}

	vector <CCAEAInd>  initial_rest_indi_pop;
	for (int n = 0; n< popsize; n++)//计算角度 第一轮按角度填充个体 并收集剩余个体
	{
		//initial_indi_pop[n].yobj2angle(ObserverPoint);
		initial_indi_pop[n].obj2angle_index(IdealPoint, popsize);

		if (sectorpopFlag[initial_indi_pop[n].sectorialindex] == false)//if (sectorialgrid[initial_indi_pop[n].sectorialindex] == -1)
		{
			population[initial_indi_pop[n].sectorialindex] = initial_indi_pop[n];
			sectorpopFlag[initial_indi_pop[n].sectorialindex] = true;
		}
		else{
			bool bReplaced;
			rest_ind = population[initial_indi_pop[n].sectorialindex];//.indiv;
			bReplaced = Pareto_HyperVolume_compare_sectorialgrid(initial_indi_pop[n]);
			if (bReplaced)
			{
				initial_rest_indi_pop.push_back(rest_ind);
			}
			else{
				initial_rest_indi_pop.push_back(initial_indi_pop[n]);
			}
		}	
	}
	initial_indi_pop.clear();

	for (int n = 0; n< popsize; n++)//填充剩余的个体
	{
		if (sectorpopFlag[n] == false)
		{
			double min_angle_difference, angle_difference;
			int min_angle_difference_ind;
			int rest_size = initial_rest_indi_pop.size();

			min_angle_difference = fabs(initial_rest_indi_pop[0].sectorialangle - (double)n / (double)(popsize - 1));
			//min_angle_difference = abs( initial_rest_indi_pop[0].sectorialindex - n );
			min_angle_difference_ind = 0;
			for (int j = 1; j<rest_size; j++)
			{
				angle_difference = fabs(initial_rest_indi_pop[j].sectorialangle - (double)n / (double)(popsize - 1));
				//angle_difference = abs( initial_rest_indi_pop[j].sectorialindex - n );
				if (angle_difference < min_angle_difference)
				{
					min_angle_difference = angle_difference;
					min_angle_difference_ind = j;
				}
			}
			population[n] = initial_rest_indi_pop[min_angle_difference_ind];
			sectorpopFlag[n] = true;
			rest_ind = initial_rest_indi_pop[min_angle_difference_ind];
			initial_rest_indi_pop[min_angle_difference_ind] = initial_rest_indi_pop[rest_size - 1];
			initial_rest_indi_pop[rest_size - 1] = rest_ind;
			initial_rest_indi_pop.pop_back();
		}
	}

	initial_rest_indi_pop.clear();
	sectorpopFlag.clear();
}

void CCAEA::execute(int mg, int irun, vector<double>& igd, vector<double>& hv, double &runtime, PTYPE ptype,MGraph *G)
{
	////seed = (seed + 23)%1377;
	////rnd_uni_init = -(long)seed;
	//vector<double> gd;
	initTopologies();
	if (numprocs >1) {
		//分配缓冲空间
		//更新理想点使用的缓冲区
		if (ptype == IDEAL_INDIV || ptype == IDEAL)
		{
			recv_idealpoint = new double*[partition_niche];
			send_idealpoint = new double*[partition_niche];
			req_idealpoint_recv = new MPI_Request[partition_niche];
			req_idealpoint_send = new MPI_Request[partition_niche];
			isBegin_Ideal = new bool[partition_niche];
		}
		//个体更新接收
		//决策向量+目标向量+权重向量
		s_indiv = 1 + nvar + nobj;
		indivBufSize = 1 + s_indiv * popsize;//sizeof(int)+(sizeof(int) + sizeof(double)*(nvar+nobj)) * maxNiche;
		if (ptype == IDEAL_INDIV || ptype == INDIV)
		{
			recv_indiv = new double*[partition_niche];
			send_indiv = new double*[partition_niche];
			req_indiv_recv = new MPI_Request[partition_niche];
			req_indiv_send = new MPI_Request[partition_niche];
			isBegin_Indiv = new bool[partition_niche];
			indiv_buf = new double[indivBufSize];
		}
		//cout << "[" << mpi_rank << "]  3" << endl;
		//cout << "[" << mpi_rank << "] partition niche :" << partition_niche << endl;
		for (int i = 0; i < partition_niche; i++)
		{
			if (ptype == IDEAL_INDIV || ptype == IDEAL)
			{
				recv_idealpoint[i] = new double[nobj];
				send_idealpoint[i] = new double[nobj];
				isBegin_Ideal[i] = false;
				//非阻塞接收理想点更新
				MPI_Irecv(recv_idealpoint[i], nobj, MPI_DOUBLE, neighbor_partitions[i], TAG_IDEALPOINT, MPI_COMM_WORLD, &(req_idealpoint_recv[i]));
			}

			if (ptype == IDEAL_INDIV || ptype == INDIV)
			{
				recv_indiv[i] = new double[indivBufSize];
				send_indiv[i] = new double[indivBufSize];
				isBegin_Indiv[i] = false;
				MPI_Irecv(recv_indiv[i], indivBufSize, MPI_DOUBLE, neighbor_partitions[i], TAG_INDIV, MPI_COMM_WORLD, &(req_indiv_recv[i]));
			}
		}
		//cout<<"["<<mpi_rank<<"] "<<endl;
		MPI_Barrier(MPI_COMM_WORLD);
	}

	isIdealpointUpdate = false;
	recordUpdateofIndiv = new bool[population.size()];
	for (int i = 0; i < population.size(); i++)
		recordUpdateofIndiv[i] = false;

	init_population(G);

	if (ptype == IDEAL_INDIV || ptype == IDEAL)
	{
		if (numprocs > 1) {
			//cout << "[" <<mpi_rank << "] sync idealpoint" <<endl;
			for (int i = 0; i < nobj; i++)
			{
				send_idealpoint[0][i] = IdealPoint[i];
			}

			//cout << "[" << mpi_rank << "] ideal (" << idealpoint[0] << "," << idealpoint[1] << "," << idealpoint[2] << ")" << endl;

			double *recv_utopian;
			if (mpi_rank == 0)
				recv_utopian = new double[nobj*numprocs];

			//rank = 0 收集理想点
			//收集数据
			int *rcounts = new int[numprocs];
			for (int i = 0; i < numprocs; i++) {
				rcounts[i] = nobj;
			}
			int *displs = new int[numprocs];
			displs[0] = 0;
			for (int i = 1; i < numprocs; i++)
				displs[i] = displs[i - 1] + rcounts[i - 1];

			MPI_Gatherv(send_idealpoint[0], nobj, MPI_DOUBLE, recv_utopian, rcounts, displs, MPI_DOUBLE, root_rank, MPI_COMM_WORLD);

			//统一理想点
			if (mpi_rank == root_rank)
			{
				for (int i = 0; i < numprocs; i++)
				{
					int j;
					for (j = 0; j < nobj; j++)
					{
						double tmp = recv_utopian[i*nobj + j];
						if (tmp < IdealPoint[j])
							IdealPoint[j] = tmp;
					}
				}
				//cout << "tongyi(" << idealpoint[0] << "," << idealpoint[1] << "," << idealpoint[2] << ")" << endl;
			}
			//分发、同步理想点
			for (int i = 0; i < nobj; i++)
			{
				send_idealpoint[0][i] = IdealPoint[i];
			}
			MPI_Bcast(send_idealpoint[0], nobj, MPI_DOUBLE, root_rank, MPI_COMM_WORLD);

			//更新进程内理想点
			for (int i = 0; i < nobj; i++) {
				IdealPoint[i] = send_idealpoint[0][i];
			}
			//cout << "[" << mpi_rank << "] after sync idealpoint("<<idealpoint[0]<<","<<idealpoint[1]<<","<<idealpoint[2]<<")" << endl;
			if (mpi_rank == root_rank)
				delete[]recv_utopian;
			delete[]rcounts;
			delete[]displs;
		}
	}


	int gen = 1;
	double igd_gen = 0;
	double hv_gen = 0;
	char filename[1024];
	clock_t begintime;
	clock_t endtime;
	unevolvetime = 0;

	if (IS_DEBUG) {
		begintime = clock();
		//收集计算解集质量的必要数据,目标向量
		//等待，同步
		if (mpi_rank == 0)
		{
			//载入真实前沿
			sprintf(filename, "%s/TestData/PF_Real/%s(%d).dat", WORKINGDIR, strTestInstance, nobj);
			loadpfront(filename, ps);
		}

		if (MODE_DEBUG)
		{
			//cout << "[" << mpi_rank << "] before gather pop y" << endl;
			gather_pop_y();
			//cout << "[" << mpi_rank << "] after" << endl;
			//Root进程处理数据
			if (mpi_rank == root_rank)
			{
				//计算解集质量（IGD，HV)
				calc_qualities(igd_gen, hv_gen);
				igd.push_back(gen);  igd.push_back(igd_gen);
				hv.push_back(gen);  hv.push_back(hv_gen);
				cout << "gen =  " << gen << "  hypervolume = " << hv_gen << "  " << "  IGD = " << igd_gen << "  " << endl;
			}
			//MPI_Barrier(MPI_COMM_WORLD);
		}
		//等待，同步
		MPI_Barrier(MPI_COMM_WORLD);
		endtime = clock();
		unevolvetime += endtime - begintime;
	}

	int iLen = 0;
	int Len = ceil(popsize / 10.0);
	CCAEAInd child2;
	//开始执行进化
	for (gen = 2; gen <= mg; gen++)
	{
		for (int i = 0; i < pops_sub_self[mpi_rank]; i++)
		{
			int parent_index1, parent_index2;

			iLen = iLen % Len;
			if (iLen == 0)
				parent_index1 = 0;
			else if (iLen == 1)
				parent_index1 = popsize - 1;
			else
				parent_index1 = tour_selection_hv2(population);// int(rnd_uni(&rnd_uni_init)*( popsize )); ////int(rnd_uni(&rnd_uni_init)*( popsize )); //
			//ind_selected = &(population[parent_index1]);

			parent_index2 = tour_selection_hv2(population);// int(rnd_uni(&rnd_uni_init)*( popsize ));// //int(rnd_uni(&rnd_uni_init)*( popsize )); //
			//ind_selected2 = &(population[parent_index2]);

			//int parent_index3 =  tour_selection_hv2( population );

			//sbx_mutation_evolution(*ind_selected, *ind_selected2);
			//if(rnd_uni(&rnd_uni_init)<0.9)
			real_sbx_xoverA(population[parent_index1], population[parent_index2], onechild, child2);
			//diff_evo_xoverB(population[parent_index1], population[parent_index2], population[parent_index3], onechild, 0.5, 1.0);
			realmutation(onechild, 1.0 / nvar);

			onechild.obj_eval();

			//nfes++;
			iLen++;

			bool updated = update_extreme_point(onechild);//update_anglerange_maxminindi(onechild, true); //update_nondominated_obj_maxmin(onechild);//update_angle_maxmin(onechild);//
			if (updated)
			{
				reset_angle();
			}

			//onechild.yobj2angle(ObserverPoint);
			onechild.obj2angle_index(IdealPoint, popsize);//之前之所以出现奇异点 是因为计算onechild角度是在统计极端点前 极端点改变了 而onechild的角度是老的

			//update_sectorialgrid(onechild);
			if (Pareto_HyperVolume_compare_sectorialgrid(onechild))
				recordUpdateofIndiv[onechild.sectorialindex] = true;
		}



		if (ptype == IDEAL_INDIV || ptype == IDEAL)
		{
			if (numprocs > 1)
			{
				cout << gen << "    [" << mpi_rank << "] update idealpoint" << endl;
				if (isIdealpointUpdate)
				{
					//cout << "has update" << endl;
					for (int i = 0; i < partition_niche; i++)
					{
						sendIdealpoint(send_idealpoint[i], isBegin_Ideal[i], neighbor_partitions[i], &(req_idealpoint_send[i]));
					}
					isIdealpointUpdate = false;
				}
				//cout << "[" << mpi_rank << "] update idealpoint finish" << endl;
			}
			update_idealpoint_from_outside();
		}

		if (ptype == IDEAL_INDIV || ptype == INDIV)
		{
			if (numprocs > 1 && gen % comm_i == 0)
			{
				cout << "[" << mpi_rank << "] update indiv" << endl;
				//更新个体

				int updated_size = 0;
				int m = 1;
				for (int i = 0; i < population.size(); i++)
				{
					if ((i < pops_sub_start[mpi_rank] || i >(pops_sub_start[mpi_rank] + pops_sub_self[mpi_rank])) && (recordUpdateofIndiv[i] == true))
					{
						indiv_buf[m++] = i;
						for (int v = 0; v < nvar; v++)
							indiv_buf[m++] = population[i].x_var[v];
						for (int v = 0; v < nobj; v++)
							indiv_buf[m++] = population[i].y_obj[v];
						updated_size++;
					}
				}

				if (updated_size > 0)
				{
					indiv_buf[0] = updated_size;
					int send_size = 1 + updated_size * s_indiv;
					for (int i = 0; i < partition_niche; i++)
					{
						int neiProc = neighbor_partitions[i];

						//memcpy(send_indiv[i],indiv_buf, send_size * sizeof(double));
						for (int k = 0; k < 1 + updated_size*(s_indiv); k++)
							send_indiv[i][k] = indiv_buf[k];
						/*
						cout << "[" << mpi_rank << "] up size:" << send_indiv[i][0] << "     ";
						for (int k = 1; k < send_indiv[i][0]*(1 + nvar + nobj); k++)
						cout << send_indiv[i][k] << ",";
						cout << endl;
						*/
						if (isBegin_Indiv[i])
						{
							//cout << gen <<"   [" << mpi_rank << "] wait indiv" << endl;
							MPI_Wait(&req_indiv_send[i], &status_indiv);
							//MPI_Cancel(&req_indiv_send[i]);
							//cout << gen<<"   [" << mpi_rank << "] wait indiv finish" << endl;
						}
						//cout << gen<< "   [" << mpi_rank << "] -->" << neiProc << "   send indiv" << endl;
						MPI_Isend(send_indiv[i], send_size, MPI_DOUBLE, neiProc, TAG_INDIV, MPI_COMM_WORLD, &(req_indiv_send[i]));
						isBegin_Indiv[i] = true;
					}
				}

				//cout<<"["<<mpi_rank<<"] 5"<<endl;

				//重置更新标志
				for (int i = 0; i < population.size(); i++)
					recordUpdateofIndiv[i] = false;

			}
			update_indiv_from_outside();
		}

		if (IS_DEBUG && gen % I_OUT == 0) {

			//获取计算解集质量的必要数据，目标向量
			//等待，同步
			//update_problem_from_outside();
			//MPI_Barrier(MPI_COMM_WORLD);
			//update_problem_from_outside();
			if (MODE_DEBUG)
			{
				//等待，同步
				//MPI_Barrier(MPI_COMM_WORLD);
				begintime = clock();
				cout << gen << "   [" << mpi_rank << "] before gatherv" << endl;
				gather_pop_y();
				cout << "[" << mpi_rank << "] after gather" << endl;
				//Root进程处理数据
				if (mpi_rank == root_rank)
				{
					//计算解集质量（IGD，HV)
					calc_qualities(igd_gen, hv_gen);
					igd.push_back(gen);  igd.push_back(igd_gen);
					hv.push_back(gen);  hv.push_back(hv_gen);
					cout << "gen =  " << gen << "  hypervolume = " << hv_gen << "  " << "  IGD = " << igd_gen << "  " << endl;
				}
				//等待，同步
				MPI_Barrier(MPI_COMM_WORLD);
				endtime = clock();
				unevolvetime += endtime - begintime;
			}

		}

	}

	if (numprocs > 1)
	{
		if (ptype == IDEAL_INDIV || ptype == INDIV)
		{
			isDoneRecv = new int[partition_niche];
			isDoneSend = new int[partition_niche];
			b_isdone = new int[partition_niche];
			req_done_recv = new MPI_Request[partition_niche];
			req_done_send = new MPI_Request[partition_niche];

			//清理操作，完成操作
			for (int i = 0; i < partition_niche; i++)
			{
				isDoneRecv[i] = 0;
				int flag;
				MPI_Irecv(&b_isdone[i], 1, MPI_INT, neighbor_partitions[i], TAG_DONE, MPI_COMM_WORLD, &req_done_recv[i]);
				if (isBegin_Indiv[i])
				{
					MPI_Test(&req_indiv_send[i], &flag, &status_indiv);
					if (flag)
						isBegin_Indiv[i] = false;
				}
				if (isBegin_Indiv[i])
					isDoneSend[i] = 0;
				else
				{
					isDoneSend[i] = 1;
					MPI_Isend(&isDoneSend[i], 1, MPI_INT, neighbor_partitions[i], TAG_DONE, MPI_COMM_WORLD, &req_done_send[i]);
				}
			}


			while (1)
			{
				for (int i = 0; i < partition_niche; i++)
				{
					int flag;
					if (!isDoneSend[i])
					{
						if (isBegin_Indiv[i])
						{
							MPI_Test(&req_indiv_send[i], &flag, &status_indiv);
							if (flag)
								isBegin_Indiv[i] = false;
						}
						if (!isBegin_Indiv[i])
						{
							isDoneSend[i] = 1;
							MPI_Isend(&isDoneSend[i], 1, MPI_INT, neighbor_partitions[i], TAG_DONE, MPI_COMM_WORLD, &req_done_send[i]);
						}
					}
					if (!isDoneRecv[i])
					{
						MPI_Test(&req_done_recv[i], &flag, &status_done);
						if (flag && b_isdone[i])
							isDoneRecv[i] = 1;
						else
						{
							if (ptype == IDEAL_INDIV)
								update_idealpoint_from_outside();

							getAndUpdateIndivs(recv_indiv[i], indivBufSize, neighbor_partitions[i], &req_indiv_recv[i]);
						}
					}

				}
				bool isFinish = true;
				for (int i = 0; i < partition_niche; i++)
				{
					if (!isDoneRecv[i])
					{
						isFinish = false;
						break;
					}
					if (!isDoneSend[i])
					{
						isFinish = false;
						break;
					}
				}
				if (isFinish)
					break;
			}
		}

		MPI_Barrier(MPI_COMM_WORLD);

		if (ptype == IDEAL_INDIV || ptype == IDEAL)
			update_idealpoint_from_outside();

		for (int i = 0; i < partition_niche; i++)
		{
			if (ptype == IDEAL_INDIV || ptype == IDEAL)
			{
				MPI_Cancel(&req_idealpoint_recv[i]);
				MPI_Request_free(&req_idealpoint_recv[i]);
				//等待，同步
				if (isBegin_Ideal[i])
				{
					int flag = 1;
					MPI_Test(&req_idealpoint_send[i], &flag, &status_idealpoint);
					if (flag == 0)
					{
						MPI_Cancel(&req_idealpoint_send[i]);
						MPI_Request_free(&req_idealpoint_send[i]);
					}
				}
			}

			if (ptype == IDEAL_INDIV || ptype == INDIV)
			{
				MPI_Cancel(&req_indiv_recv[i]);
				MPI_Request_free(&req_indiv_recv[i]);
				if (isBegin_Indiv[i])
				{
					int flag;
					MPI_Test(&req_indiv_send[i], &flag, &status_indiv);
					if (flag)
					{
						MPI_Cancel(&req_indiv_send[i]);
						MPI_Request_free(&req_indiv_send[i]);
					}
				}
			}
		}

	}

	MPI_Barrier(MPI_COMM_WORLD);
	finish = clock();
	runtime = (double)(finish - start - unevolvetime) / CLOCKS_PER_SEC;

	//cout << "[" << mpi_rank << "] idealpoint("<<IdealPoint[0]<<","<<IdealPoint[1]<<")" << endl;
	//cout << "[" << mpi_rank << "] TrueNadirpoint(" << TrueNadirPoint[0] << "," << TrueNadirPoint[1] << ")" << endl;
	//cout << "[" << mpi_rank << "] Referentpoint(" << ReferencePoint[0] << "," << ReferencePoint[1] << ")" << endl;


	if (IS_DEBUG) {
		if (!MODE_DEBUG || (MODE_DEBUG && (gen - 1) % I_OUT != 0))
		{//最后一代的数据
			//获取计算解集质量的必要数据，目标向量
			//等待，同步
			//MPI_Barrier(MPI_COMM_WORLD);
			gather_pop_y();
			//Root进程处理数据
			if (mpi_rank == 0)
			{
				//计算解集质量（IGD，HV)
				calc_qualities(igd_gen, hv_gen);
				igd.push_back(gen - 1);  igd.push_back(igd_gen);
				hv.push_back(gen - 1);  hv.push_back(hv_gen);
				//cout<<"gen =  "<<gen-1<<"  hypervolume = "<<hv_gen<<"  "<<"  IGD = "<<igd_gen<<"  "<<endl;
			}
		}
		/*
		// save the final population - X space
		//每个进程都保存，用于后面展示
		sprintf(filename, "%s/TestData/%d/POP/POP_%s_%s(%d)_%d_%dR%d.dat",WORKINGDIR,mpi_rank,STR_MOEAD,strTestInstance,nobj, popsize, total_gen, irun);
		save_front(filename);
		sprintf(filename, "%s/TestData/%d/POF/PF_%s_%s(%d)_%d_%dR%d.dat",WORKINGDIR,mpi_rank,STR_MOEAD, strTestInstance,nobj, popsize, total_gen, irun);
		vector <CIndividual>  frontpop;
		population2front(population, frontpop);
		save_population(frontpop, filename);//保存扇形存档种群的非劣前沿

		frontpop.clear();

		sprintf(filename, "%s/TestData/%d/POS/POS_%s_%s(%d)_%d_%dR%d.dat",WORKINGDIR,mpi_rank,STR_MOEAD, strTestInstance,nobj, popsize, total_gen, irun);
		save_pos(filename);
		*/
	}
	//输出结果
	//合并前沿数据
	//等待，同步
	//	MPI_Barrier(MPI_COMM_WORLD);
	//cout << "[" << mpi_rank << "] gather" << endl;
	gather_populations();
	//cout << "[" << mpi_rank << "] gather end" << endl;
	//等待，同步
	//cout<<"["<<mpi_rank<<"] send data"<<endl;
	//MPI_Barrier(MPI_COMM_WORLD);
	if (mpi_rank == 0)
	{
		//cout << "receive mark" << endl;
		//cout<<"["<<mpi_rank<<"] 处理数据"<<endl;
		//处理数据
		vector<vector<double> > pop_x;
		vector<vector<double> > pop_y;
		vector<int> pop_index;
		pop_x.reserve(popsize);
		pop_y.reserve(popsize);
		int m = 0;
		for (int i = 0; i<popsize; i++)
		{
			vector<double> x_var;
			for (int j = 0; j<nvar; j++)
				x_var.push_back(b_recv[m++]);
			pop_x.push_back(x_var);
			vector<double> y_obj;
			for (int j = 0; j<nobj; j++)
				y_obj.push_back(b_recv[m++]);
			pop_y.push_back(y_obj);
			pop_index.push_back(i);
			//cout << "(" << y_obj[0] << "," << y_obj[1] << ")" << endl;
		}

		//cout << "receive end" << endl;
		vector<int> index_split;
		for (int i = 0; i<numprocs; i++)
		{
			for (int j = 0; j<pops_sub_self[i]; j++)
				index_split.push_back(i);
		}

		//cout << "dewfrwev" << endl;
		//保存种群前沿，未清理
		sprintf(filename, "%s/TestData/POP/POP_CAEA_%s(%d)_%d_%dR%d.dat", WORKINGDIR, strTestInstance, nobj, popsize, total_gen, irun);
		save_population_front_split(pop_y, pop_index, index_split, filename);//保存扇形存档种群的非劣前沿
		pop_index.clear();
		//过滤取出前沿
		vector <vector<double> >  pf_y;
		vector<int> pf_index;
		population2front(pop_y, pf_y, pf_index);
		//保存整体前沿
		sprintf(filename, "%s/TestData/POF/PF_CAEA_%s(%d)_%d_%dR%d.dat", WORKINGDIR, strTestInstance, nobj, popsize, total_gen, irun);
		save_population_front_split(pf_y, pf_index, index_split, filename);//保存扇形存档种群的非劣前沿
		//保存决策变量
		sprintf(filename, "%s/TestData/POS/POS_CAEA_%s(%d)_%d_%dR%d.dat", WORKINGDIR, strTestInstance, nobj, popsize, total_gen, irun);
		save_pos_split(pop_x, pf_index, index_split, filename);//保存扇形存档种群的非劣前沿

		pop_y.clear();
		pf_y.clear();
		index_split.clear();
		pf_index.clear();
	}

	population.clear();

	if (numprocs > 1)
	{
		//cout <<"["<<mpi_rank<<"] finall" << endl;
		for (int i = 0; i < partition_niche; i++)
		{
			if (ptype == IDEAL_INDIV || ptype == IDEAL)
			{
				delete[]recv_idealpoint[i];
				recv_idealpoint[i] = NULL;
				delete[]send_idealpoint[i];
				send_idealpoint[i] = NULL;
			}
			if (ptype == IDEAL_INDIV || ptype == INDIV)
			{
				delete[]recv_indiv[i];
				recv_indiv[i] = NULL;
				delete[]send_indiv[i];
				send_indiv[i] = NULL;
				delete[]indiv_buf;
				indiv_buf = NULL;
			}
		}
		if (ptype == IDEAL_INDIV || ptype == IDEAL)
		{
			delete[]recv_idealpoint;
			recv_idealpoint = NULL;
			delete[]send_idealpoint;
			send_idealpoint = NULL;
			delete[]req_idealpoint_recv;
			req_idealpoint_recv = NULL;
			delete[]req_idealpoint_send;
			req_idealpoint_send = NULL;
			delete[]isBegin_Ideal;
			isBegin_Ideal = NULL;
		}
		if (ptype == IDEAL_INDIV || ptype == INDIV)
		{
			delete[]req_indiv_recv;
			req_indiv_recv = NULL;
			delete[]req_indiv_send;
			req_indiv_send = NULL;
			delete[]isBegin_Indiv;
			isBegin_Indiv = NULL;

			delete[]req_done_recv;
			req_done_recv = NULL;
			delete[]req_done_send;
			req_done_send = NULL;
			delete[]isDoneRecv;
			isDoneRecv = NULL;
			delete[]isDoneSend;
			isDoneSend = NULL;
			delete[]b_isdone;
			b_isdone = NULL;
		}

		neighbor_partitions.clear();
		//cout << "[" << mpi_rank << "]finish" << endl;
	}
}