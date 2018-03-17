#include "dmoea-high.h"
#include "../MOEAD/dmoea.h"
#include "../MOEAD/scalarfunc.h"
//------------- Parameters in MOEA/D -------------------------

TMOEAD_HIGH::TMOEAD_HIGH(int n_H, int pops, int moead_nei)
{	//cout << "[" << mpi_rank << "] construct" << endl;
	start = clock();
	this->n_H = n_H;
	popsize = pops;
	niche = moead_nei;
	maxLayer = proc_nei;
	maxNiche = 20;


	//int *rcount = new int[numprocs + 1];
	//pops_sub = new int[numprocs];
	pops_sub_self = new int[numprocs];
	//pops_sub_start = new int[numprocs];
	int mod_tmp = popsize %numprocs;
	int v_tmp = popsize / numprocs;
	//pops_sub_start[0] = 0;
	if (mod_tmp == 0) {
		//pops_sub[0] = v_tmp;
		pops_sub_self[0] = v_tmp;
	}
	else {
		//pops_sub[0] = v_tmp + 1;//(mod_tmp> 0+1 ? 1:0);
		pops_sub_self[0] = v_tmp + 1;//(mod_tmp > 0 + 1 ? 1 : 0);
	}

	//rcount[0] = 0;
	for (int i = 1; i<numprocs; i++)
	{
		//if (mod_tmp == 0) {
			//pops_sub[i] = v_tmp;
			//pops_sub_self[i] = v_tmp;
		//}
		//else {
			//pops_sub[i] = v_tmp + (mod_tmp> i ? 1 : 0);
			pops_sub_self[i] = v_tmp + (mod_tmp > i ? 1 : 0);
		//}

		//pops_sub_start[i] = pops_sub_start[i - 1] + pops_sub[i - 1];
		//rcount[i] = rcount[i - 1] + pops_sub[i - 1];
	}
	/*
	rcount[numprocs] = rcount[numprocs - 1] + pops_sub[numprocs - 1];
	pops_sub[0] += proc_nei;//(neigh_size+1)/2;
	if (pops_sub[0] > popsize)
		pops_sub[0] = popsize;
	for (int i = 1; i<numprocs - 1; i++)
	{
		pops_sub[i] += proc_nei * 2;//(neigh_size+1)/2 * 2;
		pops_sub_start[i] -= proc_nei;//(neigh_size+1)/2;
		if (pops_sub_start[i] < 0)
			pops_sub_start[i] = 0;
		if (pops_sub_start[i] + pops_sub[i] > popsize)
			pops_sub[i] = popsize - pops_sub_start[i];
	}
	pops_sub[numprocs - 1] += proc_nei;//(neigh_size+1)/2;
	pops_sub_start[numprocs - 1] -= proc_nei;//(neigh_size+1)/2;
	if (pops_sub_start[numprocs - 1] < 0)
		pops_sub_start[numprocs - 1] = 0;
	if (pops_sub_start[numprocs - 1] + pops_sub[numprocs - 1] > popsize)
		pops_sub[numprocs - 1] = popsize - pops_sub_start[numprocs - 1];



	int *pops_lsize = new int[numprocs];
	int *pops_rsize = new int[numprocs];
	for (int i = 0; i<numprocs; i++)
	{
		pops_lsize[i] = rcount[i] - pops_sub_start[i];
		pops_rsize[i] = pops_sub_start[i] + pops_sub[i] - rcount[i + 1];
	}

	LLSize = LSize = pops_lsize[mpi_rank];
	RRSize = RSize = pops_rsize[mpi_rank];

	if (mpi_rank - 1 >= 0)
		LSize += pops_rsize[mpi_rank - 1];
	if (mpi_rank + 1 < numprocs)
		RSize += pops_lsize[mpi_rank + 1];
	//cout<<"["<<mpi_rank<<"]"<<LSize<<"\t"<<RSize<<"\t"<<pops_sub_start[mpi_rank]<<"\t"<<pops_sub[mpi_rank]<<"\t"<<pops_lsize[mpi_rank]<<"\t"<<pops_rsize[mpi_rank]<<endl;
	//MPI_Barrier(MPI_COMM_WORLD);
	delete rcount;
	delete pops_lsize;
	delete pops_rsize;
	*/
	//cout<<"["<<mpi_rank<<"]"<< pops_sub_self[mpi_rank]<<"\t"<<endl;
	//MPI_Barrier(MPI_COMM_WORLD);
	int tmp = nvar + nobj;
	s_send = tmp * pops_sub_self[mpi_rank];
	s_recv = tmp * popsize;
	b_send = new double[s_send];
	if (mpi_rank == 0)
		b_recv = new double[s_recv];

	//idealpoint = new double[nobj];
	indivpoint = new CIndividual[nobj];
	// initialize ideal point
	for (int n = 0; n < nobj; n++)
	{
		idealpoint.push_back(1.0e+30);
		indivpoint[n].rnd_init();
	}
}

TMOEAD_HIGH::~TMOEAD_HIGH()
{
	//delete pops_sub;
	//delete pops_sub_start;
	delete b_send;
	if(mpi_rank == 0)
		delete b_recv;
	//idealpoint.clear();//delete [] idealpoint;
	delete [] indivpoint;
	idealpoint.clear();
	lambda.clear();
}


void TMOEAD_HIGH::init_population()
{
	//cout << "[" << mpi_rank << "] init population" << endl;
	int size = population.size();
    for(int i=0; i<size; i++)
	{
		population[i].indiv.rnd_init();
		update_idealpoint(population[i].indiv);
	}
}


// initialize a set of evely-distributed weight vectors
void TMOEAD_HIGH::init_uniformweight()
{
	/*
	int start = pops_sub_start[mpi_rank];
	int end = start + pops_sub[mpi_rank];
    for(int i=start; i< end; i++)
	{
        TSOP sop;
		sop.array.push_back(i);
		sop.array.push_back(n_H-i);
		//printf("[%d](%d,%d)\n",mpi_rank,i,sd-i);
		for(int j=0; j<sop.array.size(); j++)
			sop.namda.push_back(1.0*sop.array[j]/n_H);

		population.push_back(sop);
	}
	*/
	//产生全部权重
	//cout << "[" << mpi_rank << "] init weight";
	vector<int> coordinate(nobj, 0);
	gen_uniformweight(nobj - 1, n_H, coordinate, n_H);
	coordinate.clear();

	if (numprocs > 1)
	{
		/*
		//cout << "[" << mpi_rank << "] lambda " << lambda.size() << endl;
		//为每个进程分配权重向量
		allocateForProcess();

		//cout << "[" << mpi_rank << "] w1" << endl;

		//分配覆盖子问题权重
		partition_overlap(maxLayer);
		//cout << "[" << mpi_rank << "] w2" << endl;
		//计算各分区的邻居
		partition_neighbor();

		//cout << "[" << mpi_rank << "] w3" << endl;
		//将进程自身权重和覆盖子问题权重加入到种群中
		set<int>::iterator it = partition_indexs[mpi_rank].begin();
		set<int>::iterator it_end = partition_indexs[mpi_rank].end();
		while (it != it_end)
		{
			TSOP sop;
			sop.array = lambda[*it];
			for (int j = 0; j < sop.array.size(); j++)
				sop.namda.push_back(1.0*sop.array[j] / n_H);
			population.push_back(sop);
			map_weightIndex_popIndex.insert(pair<int, int>(*it, population.size() - 1));
			map_popIndex_weightIndex.insert(pair<int, int>(population.size() - 1, *it));
			it++;
		}
		overlap_start_index = population.size();

		//cout << "[" << mpi_rank << "] overlap start :" << overlap_start_index << endl;;
		neighbor_overlap_indexs = vector<vector<int>>(numprocs);

		set<int>::iterator it_o = partition_overlap_indexs[mpi_rank].begin();
		set<int>::iterator it_end_o = partition_overlap_indexs[mpi_rank].end();
		while (it_o != it_end_o)
		{
			TSOP sop;
			sop.array = lambda[*it_o];
			for (int j = 0; j < sop.array.size(); j++)
				sop.namda.push_back(1.0*sop.array[j] / n_H);
			population.push_back(sop);
			int index = population.size() - 1;
			map_weightIndex_popIndex.insert(pair<int, int>(*it_o, index));
			map_popIndex_weightIndex.insert(pair<int, int>(index, *it_o));

			for (int i = 1; i < sop_partitions[*it_o].size(); i++)
			{
				neighbor_overlap_indexs[sop_partitions[*it_o][0]].push_back(index);
			}
			it_o++;
		}

		maxNiche = neighbor_overlap_indexs[0].size();
		for (int i = 0; i < numprocs; i++)
		{
			if (maxNiche < neighbor_overlap_indexs[i].size())
				maxNiche = neighbor_overlap_indexs[i].size();
		}
		*/

		//线性划分
		//cout << "[" << mpi_rank << "] lambda " << lambda.size() << endl;
		//为每个进程分配权重向量
		allocateForProcess_Linear();

		//cout << "[" << mpi_rank << "] w1" << endl;

		//分配覆盖子问题权重
		partition_overlap_linear(proc_nei);
		//cout << "[" << mpi_rank << "] w2" << endl;
		//计算各分区的邻居
		partition_neighbor();

		//cout << "[" << mpi_rank << "] w3" << endl;
		//将进程自身权重和覆盖子问题权重加入到种群中
		set<int>::iterator it = partition_indexs[mpi_rank].begin();
		set<int>::iterator it_end = partition_indexs[mpi_rank].end();
		while (it != it_end)
		{
			TSOP sop;
			sop.array = lambda[*it];
			for (int j = 0; j < sop.array.size(); j++)
				sop.namda.push_back(1.0*sop.array[j] / n_H);
			population.push_back(sop);
			map_weightIndex_popIndex.insert(pair<int, int>(*it, population.size() - 1));
			map_popIndex_weightIndex.insert(pair<int, int>(population.size() - 1, *it));
			it++;
		}
		overlap_start_index = population.size();

		//cout << "[" << mpi_rank << "] overlap start :" << overlap_start_index << endl;;
		neighbor_overlap_indexs = vector<vector<int> >(numprocs);

		set<int>::iterator it_o = partition_overlap_indexs[mpi_rank].begin();
		set<int>::iterator it_end_o = partition_overlap_indexs[mpi_rank].end();
		while (it_o != it_end_o)
		{
			TSOP sop;
			sop.array = lambda[*it_o];
			for (int j = 0; j < sop.array.size(); j++)
				sop.namda.push_back(1.0*sop.array[j] / n_H);
			population.push_back(sop);
			int index = population.size() - 1;
			map_weightIndex_popIndex.insert(pair<int, int>(*it_o, index));
			map_popIndex_weightIndex.insert(pair<int, int>(index, *it_o));

			for (int i = 1; i < sop_partitions[*it_o].size(); i++)
			{
				neighbor_overlap_indexs[sop_partitions[*it_o][0]].push_back(index);
			}
			it_o++;
		}

		maxNiche = neighbor_overlap_indexs[0].size();
		for (int i = 0; i < numprocs; i++)
		{
			if (maxNiche < neighbor_overlap_indexs[i].size())
				maxNiche = neighbor_overlap_indexs[i].size();
		}
	}
	else
	{
		for (int i = 0; i < popsize; i++)
		{
			TSOP sop;
			sop.array = lambda[i];
			for (int j = 0; j < sop.array.size(); j++)
				sop.namda.push_back(1.0*sop.array[j] / n_H);
			population.push_back(sop);
		}
	}
}

//分解子问题、并计算观察向量
void TMOEAD_HIGH::gen_uniformweight(int start_obj_index, int max_value_left, vector<int> coordinate, int sd)
{
	if (0 == start_obj_index || 0 == max_value_left)
	{
		coordinate[start_obj_index] = max_value_left;
		lambda.push_back(coordinate);
		return;
	}

	for (int i = max_value_left; i >= 0; i--)
	{
		coordinate[start_obj_index] = i;
		gen_uniformweight(start_obj_index - 1, max_value_left - i, coordinate, sd);
	}
}

void TMOEAD_HIGH::allocateForProcess()
{
	vector<bool> isAllocated(popsize, false);
	sop_partitions = vector<vector<int> >(popsize);
	queue<int> inQueue;
	set<int> inSet;
	int low = n_H;
	partition_indexs = vector<set<int> >(numprocs);

	for (int i = 0; i < numprocs; i++)
	{
		int count = 0;
		while (count < pops_sub_self[i])
		{
			if (inQueue.empty())
			{
				//cout << "["<<mpi_rank<<"]empty" << endl;
				inSet.clear();
				for (int j = 0; j < popsize; j++)
					if (isAllocated[j] == false)
					{
						inQueue.push(j);
						inSet.insert(j);
						break;
					}
			}
			int index = inQueue.front();
			inQueue.pop();
			//加入映射
			if (isAllocated[index] == false)
			{
				//权重归属partition
				sop_partitions[index].push_back(i);
				partition_indexs[i].insert(index);

				isAllocated[index] = true;
				count++;
				if (low > lambda[index][nobj - 1])
					low = lambda[index][nobj - 1];
				//计算邻居加入选择队列
				vector<int> vicinities = gen_vicinities(lambda[index]);
				//cout << "size : " << vicinities.size() << endl;
				for (int j = 0; j < vicinities.size(); j++)
					if (isAllocated[vicinities[j]] == false && inSet.count(vicinities[j]) == 0)
					{
						inQueue.push(vicinities[j]);
						inSet.insert(vicinities[j]);
					}
			}
		}
		//一个进程分配完毕
		inSet.clear();
		if ((low > 0) && !inQueue.empty())
		{
			int index = inQueue.front();
			inQueue.pop();
			while (!inQueue.empty())
			{
				int index2 = inQueue.front();
				inQueue.pop();

				if (lambda[index2][0] < lambda[index][0])
					index = index2;
				else if ((lambda[index2][0] == lambda[index][0]) && (lambda[index2][nobj - 1] > lambda[index][nobj - 1]))
					index = index2;
			}
			inQueue.push(index);
			inSet.insert(index);
		}
		else
		{
			while (!inQueue.empty())
				inQueue.pop();

			for (int j = 0; j < popsize; j++)
				if (isAllocated[j] == false)
				{
					inQueue.push(j);
					inSet.insert(j);
					break;
				}
			if (low <= 0)
				low = n_H;
		}
	}
	isAllocated.clear();
	inSet.clear();
}

void TMOEAD_HIGH::allocateForProcess_Linear()
{
	sop_partitions = vector<vector<int> >(popsize);
	partition_indexs = vector<set<int> >(numprocs);

	int index = 0;
	for (int i = 0; i < numprocs; i++)
	{
		int count = 0;
		while (count < pops_sub_self[i])
		{
			//权重归属partition
			sop_partitions[index].push_back(i);
			partition_indexs[i].insert(index++);
			count++;
		}
	}
}

void TMOEAD_HIGH::partition_overlap(int maxLayer)
{
	partition_overlap_indexs = vector<set<int> >(numprocs);
	for (int i = 0; i < numprocs; i++)
	{
		set<int> lastLay(partition_indexs[i]);

		for (int l = 0; l < maxLayer;l++)
		{
			set<int> newLay;
			set<int>::iterator it = lastLay.begin();
			set<int>::iterator it_end = lastLay.end();
			while (it != it_end)
			{
				vector<int> vicinities = gen_vicinities(lambda[*it]);
				for (int k = 0; k < vicinities.size(); k++)
				{
					if (partition_indexs[i].count(vicinities[k]) == 0 && partition_overlap_indexs[i].count(vicinities[k]) == 0)
					{
						//partition_indexs[i].insert(vicinities[k]);
						partition_overlap_indexs[i].insert(vicinities[k]);
						newLay.insert(vicinities[k]);
						sop_partitions[ vicinities[k] ].push_back(i);
					}
				}
				it++;
			}
			lastLay = newLay;
			newLay.clear();
		}
		lastLay.clear();
	}
}


void TMOEAD_HIGH::partition_overlap_linear(int niche)
{
	partition_overlap_indexs = vector<set<int> >(numprocs);

	int *pops_sub_start = new int[numprocs];
	int *pops_sub_end = new int[numprocs];
	pops_sub_start[0] = 0;
	pops_sub_end[0] = pops_sub_self[0] ;
	for (int i = 1; i<numprocs; i++)
	{
		pops_sub_start[i] = pops_sub_end[i - 1];
		pops_sub_end[i] = pops_sub_end[i-1] + pops_sub_self[i];
	}

	for (int i = 0; i < numprocs; i++)
	{
		int startN = pops_sub_start[i] - niche;
		if (startN < 0) startN = 0;
		int endN = pops_sub_end[i] + niche;
		if (endN > popsize) endN = popsize;

		for (int j = startN; j < pops_sub_start[i]; j++)
		{
			partition_overlap_indexs[i].insert(j);
			sop_partitions[j].push_back(i);
		}
		for (int j = pops_sub_end[i]; j < endN; j++)
		{
			partition_overlap_indexs[i].insert(j);
			sop_partitions[j].push_back(i);
		}
	}
}

void TMOEAD_HIGH::partition_neighbor()
{
	//cout << "[" << mpi_rank << "] : " << sop_partitions.size() << endl;
	//vector<vector<int> > nicheC(numprocs,vector<int>(numprocs,0));
	for (int i = 0; i < sop_partitions.size(); i++)
	{
		//cout << "[" << mpi_rank << "] sop : " << sop_partitions[i].size() << endl;
		//sop_partitions[i].size() > 1//有共用，即是互相覆盖区域
		for (int j = 1; j < sop_partitions[i].size(); j++)
		{
			//nicheC[sop_partitions[i][0]][sop_partitions[i][j]]++;
			if (sop_partitions[i][0] == mpi_rank)
			{
				if (neighbor_partitions_set.count(sop_partitions[i][j]) == 0)
					neighbor_partitions_set.insert(sop_partitions[i][j]);
			}
		}
	}

	//maxNiche = nicheC[0][0];
//	for each (auto arr in nicheC)
	//{
	//	for each (auto c in arr)
	//	{
	//		if (c > maxNiche)
	//			maxNiche = c;
	//	}
	//}

	//cout << "[" << mpi_rank << "] nei set size : " << neighbor_partitions_set.size() << endl;
	set<int>::iterator it = neighbor_partitions_set.begin();
	set<int>::iterator it_end = neighbor_partitions_set.end();

	while (it != it_end)
	{
		neighbor_partitions.push_back(*it);
		it++;
	}
	neighbor_partitions_set.clear();
	partition_niche = neighbor_partitions.size();
	//for (int i = 0; i < partition_niche;i++)
	//	cout << "[" << mpi_rank << "]   " << neighbor_partitions[i] << endl;
}

int TMOEAD_HIGH::cal_index_from_lambda(const vector<int> &V_obj)
{
	vector<int> h(nobj - 1, 0);
	h[nobj - 2] = n_H - V_obj[nobj - 1];
	for (int i = nobj - 3; i >0; i--)
	{
		h[i] = h[i + 1] - V_obj[i + 1];
	}
	h[0] = V_obj[0];

	int result_index = 0;
	for (int i = 0; i < nobj - 1; i++)
	{
		result_index += choose(h[i] + i, i + 1);
	}
	h.clear();

	return result_index;
}

vector<int> TMOEAD_HIGH::gen_vicinities(vector<int> V_obj)
{
	vector<int> vicinities;
	for (int i = 0; i < nobj; i++)
	{
		V_obj[i] = V_obj[i] + 1;
		if (V_obj[i] <= n_H)
		{
			for (int j = i + 1; j < nobj; j++)
			{
				V_obj[j] = V_obj[j] - 1;
				if (V_obj[j] >= 0)
				{
					int index = cal_index_from_lambda(V_obj);
					vicinities.push_back(index);
				}
				V_obj[j] = V_obj[j] + 1;
			}
		}
		V_obj[i] = V_obj[i] - 2;
		if (V_obj[i] >= 0)
		{
			for (unsigned int j = i + 1; j < nobj; j++)
			{
				V_obj[j] = V_obj[j] + 1;
				if (V_obj[j] <= n_H)
				{
					int index = cal_index_from_lambda(V_obj);
					vicinities.push_back(index);
				}
				V_obj[j] = V_obj[j] - 1;
			}
		}
		V_obj[i] = V_obj[i] + 1;
	}
	return vicinities;
}


// initialize the neighborhood of subproblems based on the distances of weight vectors
void TMOEAD_HIGH::init_neighborhood()
{
	int size = population.size();
	double *x = new double[size];
	int *idx = new int[size];
	for(int i=0; i < size; i++)
	{
		for(int j=0; j < size; j++)
		{
		    x[j]    = dist_vector(population[i].namda,population[j].namda);
			idx[j]  = j;
		}
		minfastsort(x,idx,size,niche);
		for(int k=0; k<niche; k++){
			population[i].table.push_back(idx[k]);
			//printf("[%d]%d--%d\n",mpi_rank,i,idx[k]);
		}
	}
    delete [] x;
	delete [] idx;
}

// update the best solutions of neighboring subproblems
void TMOEAD_HIGH::update_problem(CIndividual &indiv, int id)
{
	for (int i = 0; i<niche; i++)
	{
		int    k = population[id].table[i];
		double f1, f2;
		f1 = scalar_func(population[k].indiv.y_obj, population[k].namda,idealpoint, indivpoint);
		f2 = scalar_func(indiv.y_obj, population[k].namda,idealpoint, indivpoint);
		if (f2<f1){
			population[k].indiv = indiv;

			recordUpdateofIndiv[k] = true;
			/*
			update_problem_from_outside();
			if(numprocs >1)
			{
				if(k < LSize)
					sendIndiv(lsend_indiv,k,prevRank,&req_indiv[L_SEND],isLSend);
				if(k >= pops_sub[mpi_rank]-RSize)
					sendIndiv(rsend_indiv,k,nextRank,&req_indiv[R_SEND],isRSend);
			}
			*/
		}
	}
}

bool TMOEAD_HIGH::update_idealpoint_from_outside()
{
	bool isUpdate = false;
	if(numprocs >1){
		//只有多进程才使用，保证单进程可以执行
		for (int i = 0; i < partition_niche; i++)
		{
			int flag = 0;
			bool tmpUpdate = false;
			MPI_Test(&req_idealpoint_recv[i], &flag, &status_idealpoint);
			while (flag != 0)
				//if(flag != 0)
			{//接收到来自其他进程的理想点更新

				//cout << "[" << mpi_rank << "]("<<recv_idealpoint[i][0]<<","<<recv_idealpoint[i][1]<<","<<recv_idealpoint[i][2]<<")" << endl;
				tmpUpdate = update_idealpoint(recv_idealpoint[i]);
				//cout<<"["<<mpi_rank<<"]("<<recv_idealpoint[0]<<","<<recv_idealpoint[1]<<") old("<<ij[0]<<","<<ij[1]<<")now("<<idealpoint[0]<<","<<idealpoint[1]<<")"<<endl;

				//继续等待接收更新
				MPI_Irecv(recv_idealpoint[i], nobj, MPI_DOUBLE, neighbor_partitions[i], TAG_IDEALPOINT, MPI_COMM_WORLD, &req_idealpoint_recv[i]);
				flag = 0;
				MPI_Test(&req_idealpoint_recv[i], &flag, &status_idealpoint);
			}
			isUpdate |= tmpUpdate;
		}

	}
	return isUpdate;
}

void TMOEAD_HIGH::sendIdealpoint(double *buf,bool &isSend,int target,MPI_Request* request)
{
	if(numprocs >1)
	{
		//多进程才使用，保证单进程也可以执行
		//如果理想点更新，则向下一个进程发送更新信息
		//阻塞，等待上一个理想点更新信息传送完成

		//可以发送更新信息
		for(int i=0;i<nobj;i++)
		{
			buf[i] = idealpoint[i];
		}
		MPI_Isend(buf,nobj,MPI_DOUBLE,target,TAG_IDEALPOINT,MPI_COMM_WORLD,request);
		//发送指令完毕，非阻塞，继续执行接下来计算,边传输边计算
		//cout<<"["<<mpi_rank<<"] send"<<endl;
		isSend = true;
	}
}


// update the reference point
bool TMOEAD_HIGH::update_idealpoint(CIndividual &ind)
{
	bool isUpdate = false;
	for(int n=0; n<nobj; n++)
	{
		if(ind.y_obj[n]<idealpoint[n])
		{
			idealpoint[n]  = ind.y_obj[n];
			indivpoint[n]  = ind;
			isUpdate = true;
		}
	}
	if (isUpdate == true)
	{
		//for (int i = 0; i < partition_niche; i++)
			isIdealpointUpdate = true;
	}
	return isUpdate;
}


// update the reference point
bool TMOEAD_HIGH::update_idealpoint(const double* y_obj)
{
	bool isUpdate = false;
	for (int n = 0; n<nobj; n++)
	{
		if (y_obj[n]<idealpoint[n])
		{
			idealpoint[n] = y_obj[n];
			//indivpoint[n] = ind;
			isUpdate = true;
		}
	}
	//isUpdate |= update_idealpoint_from_outside();
	if (isUpdate == true)
	{
		//for (int i = 0; i < partition_niche; i++)
			//isIdealpointUpdate[i] = true;
		isIdealpointUpdate = true;
	}
	return isUpdate;
}

// recombination, mutation, update in MOEA/D
void TMOEAD_HIGH::evolution()
{
	CIndividual child, child2;
	for(int i=0; i< pops_sub_self[mpi_rank]; i++)
	{
		int   n  =  i;
		int   s  = population[n].table.size();
		int   r1 = int(s*rnd_uni(&rnd_uni_init));
		int   r2 = int(s*rnd_uni(&rnd_uni_init));
		int   p1 = population[n].table[r1];
		int   p2 = population[n].table[r2];

		//realbinarycrossover(population[p1].indiv,population[p2].indiv,child, child2);
		//diff_evo_xoverB( population[i].indiv, population[p1].indiv, population[p2].indiv, child, 0.5, 1.0);
		real_sbx_xoverA(population[p1].indiv,population[p2].indiv,child, child2);
		realmutation(child, 1.0/nvar);
		child.obj_eval();
		update_idealpoint(child);
		//if(update_idealpoint(child))
			//sendIdealpoint();
		//update_problem_from_outside();
		update_problem(child, n);
	}
}

void  TMOEAD_HIGH::execute(int mg, int irun, vector<double>& igd, vector<double>&hv, double &runtime, PTYPE ptype)
{
	// mg: maximal number of generations
	init_uniformweight();
	init_neighborhood();

	//if (mpi_rank == root_rank)
	//cout<< "["<<mpi_rank<<"] maxNiche : " << maxNiche << endl;

	//cout << "[" << mpi_rank << "] neighbor size : " << partition_niche << endl;

	if (numprocs > 1) {
		//cout << "[" << mpi_rank << "] something init";
		//只有多进程才使用，保证单进程也可以执行
		//分配缓冲空间
		//更新理想点使用的缓冲区
		if (ptype == IDEAL_INDIV || ptype == IDEAL)
		{
			recv_idealpoint = new double*[partition_niche];
			send_idealpoint = new double*[partition_niche];
			req_idealpoint_recv = new MPI_Request[partition_niche];
			req_idealpoint_send = new MPI_Request[partition_niche];
			//isIdealpointUpdate = new bool[partition_niche];
			isBegin_Ideal = new bool[partition_niche];
		}
		//个体更新接收
		//权重index + 决策向量+目标向量
		s_indiv = 1 + nvar + nobj;
		indivBufSize = 1 + s_indiv * population.size();//sizeof(int)+(sizeof(int) + sizeof(double)*(nvar+nobj)) * maxNiche;
		if (ptype == IDEAL_INDIV || ptype == INDIV)
		{
			recv_indiv = new double*[partition_niche];
			send_indiv = new double*[partition_niche];
			req_indiv_recv = new MPI_Request[partition_niche];
			req_indiv_send = new MPI_Request[partition_niche];
			isBegin_Indiv = new bool[partition_niche];
		}
			//cout << "[" << mpi_rank << "] partition niche :" << partition_niche << endl;
		for (int i = 0; i < partition_niche; i++)
		{
			if (ptype == IDEAL_INDIV || ptype == IDEAL)
			{
				recv_idealpoint[i] = new double[nobj];
				send_idealpoint[i] = new double[nobj];
				isBegin_Ideal[i] = false;
				//非阻塞接收理想点更新
				MPI_Irecv(recv_idealpoint[i], nobj, MPI_DOUBLE, neighbor_partitions[i], TAG_IDEALPOINT, MPI_COMM_WORLD, &req_idealpoint_recv[i]);
			}

			if (ptype == IDEAL_INDIV || ptype == INDIV)
			{
				recv_indiv[i] = new double[indivBufSize];
				send_indiv[i] = new double[1+s_indiv*maxNiche];
				isBegin_Indiv[i] = false;
				MPI_Irecv(recv_indiv[i], indivBufSize, MPI_DOUBLE, neighbor_partitions[i], TAG_INDIV, MPI_COMM_WORLD, &req_indiv_recv[i]);
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	isIdealpointUpdate = false;
	recordUpdateofIndiv = new bool[population.size()];
	for (int i = 0; i < population.size(); i++)
		recordUpdateofIndiv[i] = false;

	//cout<<"["<<mpi_rank<<"] 4 before"<<endl;
	init_population();
	//cout<<"["<<mpi_rank<<"] 4 after"<<population.size()<<endl;

	if (ptype == IDEAL_INDIV || ptype == IDEAL)
	{
		if (numprocs > 1) {
			for (int i = 0; i < nobj; i++)
			{
				send_idealpoint[0][i] = idealpoint[i];
			}

			//cout << "[" << mpi_rank << "] ideal (" << idealpoint[0] << "," << idealpoint[1] << "," << idealpoint[2] << ")" << endl;

			double *recv_utopian;
			if (mpi_rank == 0)
				recv_utopian = new double[nobj*numprocs];

			//rank=0收集理想点
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
						if (tmp < idealpoint[j])
							idealpoint[j] = tmp;
					}
				}
				//cout << "tongyi(" << idealpoint[0] << "," << idealpoint[1] << "," << idealpoint[2] << ")" << endl;
			}
			//分发、同步理想点
			for (int i = 0; i < nobj; i++)
			{
				send_idealpoint[0][i] = idealpoint[i];
			}
			MPI_Bcast(send_idealpoint[0], nobj, MPI_DOUBLE, root_rank, MPI_COMM_WORLD);

			//更新进程内理想点
			for (int i = 0; i < nobj; i++) {
				idealpoint[i] = send_idealpoint[0][i];
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
			if (mpi_rank == 0)
			{
				//计算解集质量（IGD，HV)
				calc_qualities(igd_gen, hv_gen);
				igd.push_back(gen);  igd.push_back(igd_gen);
				hv.push_back(gen);  hv.push_back(hv_gen);
				cout << "gen =  " << gen << "  hypervolume = " << hv_gen << "  " << "  IGD = " << igd_gen << "  " << endl;
			}
		}
		//等待，同步
		MPI_Barrier(MPI_COMM_WORLD);
		endtime = clock();
		unevolvetime += endtime - begintime;
	}

	//开始执行进化
	for (gen = 2; gen <= mg; gen++)
	{
		//cout<<"["<<mpi_rank<<"] gen"<<gen<<endl;
		evolution();

		if (ptype == IDEAL_INDIV || ptype == IDEAL)
		{
			if (numprocs > 1)
			{
				//cout << "[" << mpi_rank << "] update idealpoint" << endl;
				if (isIdealpointUpdate)
				{
					//cout << "has update" << endl;
					for (int i = 0; i < partition_niche; i++)
					{
						sendIdealpoint(send_idealpoint[i], isBegin_Ideal[i], neighbor_partitions[i], &req_idealpoint_send[i]);
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
				//cout << "[" << mpi_rank << "] update indiv begin" << endl;
				//cout << "[" << mpi_rank << "] update indiv" << endl;
				//更新个体
				for (int i = 0; i < partition_niche; i++)
				{
					vector<int> updated_indivs;
					int neiProc = neighbor_partitions[i];
					for (int j = 0; j < neighbor_overlap_indexs[neiProc].size(); j++)
					{
						if (recordUpdateofIndiv[neighbor_overlap_indexs[neiProc][j]])
							updated_indivs.push_back(neighbor_overlap_indexs[neiProc][j]);
					}
					int updated_size = updated_indivs.size();

					if (updated_size > 0)
					{
						if (isBegin_Indiv[i])
						{
							//cout << "[" << mpi_rank << "] wait indiv" << endl;
							MPI_Wait(&req_indiv_send[i], &status_indiv);
							//cout << "[" << mpi_rank << "] wait indiv finish" << endl;
						}
						int m = 0;
						send_indiv[i][m++] = updated_size;
						for (int j = 0; j < updated_size; j++)
						{
							send_indiv[i][m++] = map_popIndex_weightIndex[updated_indivs[j]];
							int popIndex = updated_indivs[j];
							for (int v = 0; v < nvar; v++)
								send_indiv[i][m++] = population[popIndex].indiv.x_var[v];
							for (int v = 0; v < nobj; v++)
								send_indiv[i][m++] = population[popIndex].indiv.y_obj[v];
						}

						MPI_Isend(send_indiv[i], 1 + s_indiv*updated_size, MPI_DOUBLE, neighbor_partitions[i], TAG_INDIV, MPI_COMM_WORLD, &req_indiv_send[i]);
						isBegin_Indiv[i] = true;
					}
					updated_indivs.clear();
					//cout << "[" << mpi_rank << "]--" << neighbor_partitions[i] << "--update size : " << updated_size << endl;
				}
				//cout << "[" << mpi_rank << "] update indiv finish" << endl;
				//cout<<"["<<mpi_rank<<"] 5"<<endl;

				//重置更新标志
				for (int i = 0; i < population.size(); i++)
					recordUpdateofIndiv[i] = false;

				//cout << "[" << mpi_rank << "] update indiv end" << endl;
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
				gather_pop_y();

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
							if(ptype == IDEAL_INDIV )
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

		if(ptype == IDEAL_INDIV || ptype == IDEAL)
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

	if(IS_DEBUG){
		if(!MODE_DEBUG || (MODE_DEBUG && (gen-1)%I_OUT !=0))
		{//最后一代的数据
			//获取计算解集质量的必要数据，目标向量
			//等待，同步
			//MPI_Barrier(MPI_COMM_WORLD);
			gather_pop_y();
			//Root进程处理数据
			if(mpi_rank ==0)
			{
				//计算解集质量（IGD，HV)
				calc_qualities(igd_gen,hv_gen);
				igd.push_back(gen-1);  igd.push_back(igd_gen);
				hv.push_back(gen-1);  hv.push_back(hv_gen);
				cout<<"Final : "<<gen-1<<"  hypervolume = "<<hv_gen<<"  "<<"  IGD = "<<igd_gen<<"  "<<endl;
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
	if(mpi_rank == 0)
	{
		//cout << "receive mark" << endl;
		//cout<<"["<<mpi_rank<<"] 处理数据"<<endl;
		//处理数据
		vector<vector<double> > pop_x;
		vector<vector<double> > pop_y;
		vector<int> pop_index;
		pop_x.reserve(popsize);
		pop_y.reserve(popsize);
		int m=0;
		for(int i=0;i<popsize;i++)
		{
			vector<double> x_var;
			for(int j=0;j<nvar;j++)
				x_var.push_back(b_recv[m++]);
			pop_x.push_back(x_var);
			vector<double> y_obj;
			for(int j=0;j<nobj;j++)
				y_obj.push_back(b_recv[m++]);
			pop_y.push_back(y_obj);
			pop_index.push_back(i);
			//cout << "(" << y_obj[0] << "," << y_obj[1] << ")" << endl;
		}

		//cout << "receive end" << endl;
		vector<int> index_split;
		for(int i =0;i<numprocs;i++)
		{
			for(int j = 0;j<pops_sub_self[i];j++)
				index_split.push_back(i);
		}

		//cout << "dewfrwev" << endl;
		//保存种群前沿，未清理
		sprintf(filename, "%s/TestData/POP/POP_%s_%s(%d)_%d_%dR%d.dat",WORKINGDIR,STR_MOEAD, strTestInstance,nobj, popsize, total_gen, irun);
		save_population_front_split(pop_y,pop_index,index_split,filename);//保存扇形存档种群的非劣前沿
		pop_index.clear();
		//过滤取出前沿
		vector <vector<double> >  pf_y;
		vector<int> pf_index;
		population2front(pop_y,pf_y,pf_index);
		//保存整体前沿
		sprintf(filename, "%s/TestData/POF/PF_%s_%s(%d)_%d_%dR%d.dat",WORKINGDIR,STR_MOEAD, strTestInstance,nobj, popsize, total_gen, irun);
		save_population_front_split(pf_y,pf_index,index_split,filename);//保存扇形存档种群的非劣前沿
		//保存决策变量
		sprintf(filename, "%s/TestData/POS/POS_%s_%s(%d)_%d_%dR%d.dat",WORKINGDIR,STR_MOEAD, strTestInstance,nobj, popsize, total_gen, irun);
		save_pos_split(pop_x,pf_index,index_split,filename);//保存扇形存档种群的非劣前沿

		pop_y.clear();
		pf_y.clear();
		index_split.clear();
		pf_index.clear();
	}

	population.clear();

	if(numprocs > 1)
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
		sop_partitions.clear();
		partition_indexs.clear();
		partition_overlap_indexs.clear();
		neighbor_partitions_set.clear();
		neighbor_partitions.clear();
		map_weightIndex_popIndex.clear();
		map_popIndex_weightIndex.clear();
		neighbor_overlap_indexs.clear();
		//cout << "[" << mpi_rank << "]finish" << endl;
	}
}

void TMOEAD_HIGH::save_front(char saveFilename[1024])
{
    std::fstream fout;
	fout.open(saveFilename,std::ios::out);
	for(int n=0; n<population.size(); n++)
	{
		for(int k=0;k<nobj;k++)
			fout<<population[n].indiv.y_obj[k]<<"  ";
		fout<<"\n";
	}
	fout.close();
}

void TMOEAD_HIGH::save_pos(char saveFilename[1024])
{
	std::fstream fout;
	fout.open(saveFilename,std::ios::out);
	for(int n=0; n<population.size(); n++)
	{
		for(int k=0;k<nvar;k++)
			fout<<population[n].indiv.x_var[k]<<"  ";
		fout<<"\n";
	}
	fout.close();
}


void TMOEAD_HIGH::operator=(const TMOEAD_HIGH &emo)
{
    popsize        = emo.popsize;
	population  = emo.population;
	indivpoint  = emo.indivpoint;
	niche       = emo.niche;
}

void TMOEAD_HIGH::gather_pop_y()
{
	//发送各个个体的目标向量，用于计算IGD和HV指标
	int m=0;
	for(int i =0;i<pops_sub_self[mpi_rank];i++)
	{
	 //for(int j =0;j<nvar;j++)
		//b_send[m++] = population[i].indiv.x_var[j];
		for(int j=0;j<nobj;j++){
			b_send[m++] = population[i].indiv.y_obj[j];
		}
	}
	int size_send = nobj * pops_sub_self[mpi_rank];

	//cout << "["<<mpi_rank<<"]  "<<size_send<<"\t"<<m << endl;
	//收集数据
	int *rcounts = new int[numprocs];
	for(int i=0;i<numprocs;i++)
		rcounts[i] = nobj*pops_sub_self[i];
	int *displs = new int[numprocs];
	displs[0] = 0;
	for(int i=1;i<numprocs;i++)
		displs[i] = displs[i-1] + rcounts[i-1];

	//cout << "[" << mpi_rank << "] before gatherv" << endl;
	MPI_Gatherv(b_send,size_send,MPI_DOUBLE,b_recv,rcounts,displs,MPI_DOUBLE, root_rank,MPI_COMM_WORLD);
	delete rcounts;
	rcounts = NULL;
	delete displs;
	displs = NULL;
}

void TMOEAD_HIGH::gather_populations()
{
	//cout << "["<<mpi_rank<<"]mark" << endl;
	//发送各个个体的目标向量，用于计算IGD和HV指标
	int m = 0;
	for (int i = 0; i < pops_sub_self[mpi_rank]; i++)
	{
		for (int j = 0; j < nvar; j++)
			b_send[m++] = population[i].indiv.x_var[j];
		for (int j = 0; j < nobj; j++) {
			b_send[m++] = population[i].indiv.y_obj[j];
		}
	}
	//cout <<"["<<mpi_rank<<"]"<<s_send<< "  adasd:"<<LLSize << "\t" << LLSize + pops_sub_self[mpi_rank] << endl;
	//收集数据
	int *rcounts = new int[numprocs];
	for(int i=0;i<numprocs;i++)
		rcounts[i] = (nvar+nobj)*pops_sub_self[i];
	int *displs = new int[numprocs];
	displs[0] = 0;
	for(int i=1;i<numprocs;i++)
		displs[i] = displs[i-1] + rcounts[i-1];
	MPI_Gatherv(b_send,s_send,MPI_DOUBLE,b_recv,rcounts,displs,MPI_DOUBLE, root_rank,MPI_COMM_WORLD);
	delete rcounts;
	rcounts = NULL;
	delete displs;
	displs = NULL;
	//cout << "["<<mpi_rank<<"]mark end" << endl;
}

void TMOEAD_HIGH::calc_qualities(double &igd,double& hv)
{//只有Root进程才使用
	//处理数据
	vector<vector<double> > pop_y;
	pop_y.reserve(popsize);
	int m=0;
	for(int i=0;i<popsize;i++)
	{
		vector<double> y_obj;
		for(int j=0;j<nobj;j++)
			y_obj.push_back(b_recv[m++]);
		//printf("[%d](%f,%f)\n",mpi_rank,y_obj[0],y_obj[1]);
		pop_y.push_back(y_obj);
	}

	//过滤取出前沿
	vector <vector<double> >  pf_y;
	population2front(pop_y,pf_y);
	//计算IGD
	igd = calc_distance(pf_y);
	//计算HV
	hv = calc_hypervolume(pf_y);
	pop_y.clear();
	pf_y.clear();
}

double TMOEAD_HIGH::calc_distance(vector<vector<double> >& ParetoFront)
{
	double distance = 0;
	for(int i=0; i<ps.size(); i++)
	{
		double min_d = 1.0e+10;
		for(int j=0; j<ParetoFront.size(); j++)
		{
			double d = dist_vector(ps[i].y_obj, ParetoFront[j]);
			if(d<min_d)  min_d = d;
		}
		distance+= min_d;
	}
	distance/=ps.size();

	return distance;
}

double TMOEAD_HIGH::calc_hypervolume(vector<vector<double> >& ParetoFront)
{
	wfg hv;
	double result = hv.compute(ParetoFront,hv_ref_ub);

	return result / cubeHV;
}

void TMOEAD_HIGH::getAndUpdateIndivs(double* b_recv_indiv,int datasize,int source,MPI_Request *request)
{
	int flag =0;
	//cout << "[" << mpi_rank << "] update indiv start" << endl;
	MPI_Test(request,&flag,&status_indiv);
	if(flag !=0)
	{//接收到来自进程的个体更新
		//从缓存中获取数据
		//cout << "[" << mpi_rank << "] receive" << endl;
		int m = 0;
		int updated_size =(int) b_recv_indiv[m++];
		for (int i = 0; i < updated_size; i++)
		{
			CIndividual indiv;
			int weightIndex = (int)b_recv_indiv[m++];
			for (int v = 0; v < nvar; v++)
				indiv.x_var[v] = b_recv_indiv[m++];
			for (int v = 0; v < nobj; v++)
				indiv.y_obj[v] = b_recv_indiv[m++];

			//cout << "[" << mpi_rank << "] " << weightIndex <<"("<<lambda[weightIndex][0]<<","<<lambda[weightIndex][1]<<","<<lambda[weightIndex][2]<<")"<< endl;
			//cout << "[" << mpi_rank << "] " << indiv.x_var[0] << "," << indiv.x_var[1] << "," << indiv.x_var[2] << "," << indiv.x_var[3] << endl;
			//cout << "[" << mpi_rank << "](" << indiv.y_obj[0] << ","<<indiv.y_obj[1]<<","<<indiv.y_obj[2] << ")" << endl;
			//更新理想点
			update_idealpoint(indiv);

			//更新个体
			int index = map_weightIndex_popIndex[weightIndex];
			//cout << "[" << mpi_rank << "] pop : " << index << endl;
			for (int j = 0; j < niche; j++)
			{
				int k = population[index].table[j];
				double f1, f2;
				f1 = scalar_func(population[k].indiv.y_obj, population[k].namda, idealpoint, indivpoint);
				f2 = scalar_func(indiv.y_obj, population[k].namda, idealpoint, indivpoint);
				if (f2 < f1) {
					population[k].indiv = indiv;
					recordUpdateofIndiv[k] = false;
				}
			}
		}
		//cout << "[" << mpi_rank << "]  asda" << endl;
		MPI_Irecv(b_recv_indiv, datasize, MPI_DOUBLE, source, TAG_INDIV, MPI_COMM_WORLD, request);
		//cout << "[" << mpi_rank << "] after" << endl;
	}
}

void TMOEAD_HIGH::sendIndiv(double* b_send_indiv,int index,int target,MPI_Request *request,bool &isBegin)
{
	//等待上一个理想点更新信息传送完成
	int flag = 1;
	if(isBegin)
	{
		//MPI_Test(request,&flag,&status_indiv);
		//if(flag == 0)
		{
			//因为个体一定比以前的好
			//MPI_Cancel(request);
		}
		MPI_Wait(request,&status_indiv);
	}

	//可以发送更新信息
	int m=0;
	TSOP &targetP = population[index];
	for(int i=0;i<nvar;i++)
		b_send_indiv[m++] = targetP.indiv.x_var[i];
	for(int i=0;i<nobj;i++)
		b_send_indiv[m++] = targetP.indiv.y_obj[i];
	for(int i=0;i<nobj;i++)
		b_send_indiv[m++] = targetP.namda[i];

	MPI_Isend(b_send_indiv,s_indiv,MPI_DOUBLE,target,TAG_INDIV,MPI_COMM_WORLD,request);
	//发送指令完毕，非阻塞，继续执行接下来计算,边传输边计算
	//cout<<"["<<mpi_rank<<"] send"<<endl;
	isBegin = true;
}

void TMOEAD_HIGH::update_indiv_from_outside()
{
	if (numprocs > 1)
	{
		for (int i = 0; i < partition_niche;i++)
			getAndUpdateIndivs(recv_indiv[i],indivBufSize,neighbor_partitions[i], &req_indiv_recv[i]);
	}
}
