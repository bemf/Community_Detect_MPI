#include "dmoea.h"
#include "scalarfunc.h"
//------------- Parameters in MOEA/D -------------------------

TMOEAD::TMOEAD(int n_H,int pops,int moead_nei)
{
	start = clock();
	this->n_H = n_H;
	popsize = pops;
	niche = moead_nei;

	int *rcount = new int[numprocs + 1];
	pops_sub = new int[numprocs];
	pops_sub_self = new int[numprocs];
	pops_sub_start = new int[numprocs];
	int mod_tmp = popsize %numprocs;
	int v_tmp = popsize / numprocs;
	pops_sub_start[0] = 0;
	if (mod_tmp == 0) {
		pops_sub[0] = v_tmp;
		pops_sub_self[0] = v_tmp;
	}
	else {
		pops_sub[0] = v_tmp + 1;//(mod_tmp> 0+1 ? 1:0);
		pops_sub_self[0] = v_tmp + 1;//(mod_tmp > 0 + 1 ? 1 : 0);
	}
	rcount[0] = 0;
	for (int i = 1; i<numprocs; i++)
	{
		if (mod_tmp == 0) {
			pops_sub[i] = v_tmp;
			pops_sub_self[i] = v_tmp;
		}
		else {
			pops_sub[i] = v_tmp + (mod_tmp> i ? 1 : 0);
			pops_sub_self[i] = v_tmp + (mod_tmp > i ? 1 : 0);
		}

		pops_sub_start[i] = pops_sub_start[i - 1] + pops_sub[i - 1];
		rcount[i] = rcount[i - 1] + pops_sub[i - 1];
	}
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

	int tmp = nvar+nobj;
	s_send = tmp * pops_sub_self[mpi_rank];
	s_recv = tmp * popsize;
	b_send = new double[s_send];
	if(mpi_rank == 0)
		b_recv = new double[s_recv];

	//idealpoint = new double[nobj];
	indivpoint = new CIndividual[nobj];
	// initialize ideal point
    for(int n=0; n<nobj; n++)
	{
		idealpoint.push_back( 1.0e+30 );
		indivpoint[n].rnd_init();
	}

	isIdealpointUpdate[LEFT] = false;
	isIdealpointUpdate[RIGHT] = false;

	recordUpdateofIndiv = new bool[ pops_sub[mpi_rank] ];
	//重置更新标志
	for(int i=0;i<pops_sub[mpi_rank];i++)
		recordUpdateofIndiv[i] = false;
}

TMOEAD::~TMOEAD()
{
	delete pops_sub;
	delete pops_sub_start;
	delete b_send;
	if(mpi_rank == 0)
		delete b_recv;
	//idealpoint.clear();//delete [] idealpoint;
	delete [] indivpoint;
	idealpoint.clear();

}


void TMOEAD::init_population()
{
    for(int i=0; i<pops_sub[mpi_rank]; i++)
	{
		population[i].indiv.rnd_init();
		update_idealpoint(population[i].indiv);
		//if(update_idealpoint(population[i].indiv))
			//sendIdealpoint();
	}
}


// initialize a set of evely-distributed weight vectors
void TMOEAD::init_uniformweight()
{
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
}

//分解子问题、并计算观察向量
void TMOEAD::gen_uniformweight(int start_obj_index, int max_value_left, vector<int> coordinate, int sd)
{
	if (0 == start_obj_index || 0 == max_value_left)
	{
		coordinate[start_obj_index] = max_value_left;
		TSOP sop;
		for (int k = 0; k < nobj; k++)
		{
			sop.array.push_back(coordinate[k]);
			sop.namda.push_back(1.0*sop.array[k] / sd);
		}
		population.push_back(sop);
		return;
	}

	for (int i = max_value_left; i >= 0; i--)
	{
		coordinate[start_obj_index] = i;
		gen_uniformweight(start_obj_index - 1, max_value_left - i, coordinate, sd);
	}
}

// initialize the neighborhood of subproblems based on the distances of weight vectors
void TMOEAD::init_neighborhood()
{
    double *x   = new double[pops_sub[mpi_rank]];
	int    *idx = new int[pops_sub[mpi_rank]];
	niche = niche <= pops_sub[mpi_rank] ? niche : pops_sub[mpi_rank];
	for(int i=0; i<pops_sub[mpi_rank]; i++)
	{
		for(int j=0; j<pops_sub[mpi_rank]; j++)
		{
		    x[j]    = dist_vector(population[i].namda,population[j].namda);
			idx[j]  = j;
		}
		minfastsort(x,idx,pops_sub[mpi_rank],niche);
		for(int k=0; k<niche; k++){
			population[i].table.push_back(idx[k]);
			//printf("[%d]%d--%d\n",mpi_rank,i,idx[k]);
		}
	}
    delete [] x;
	delete [] idx;
}

// update the best solutions of neighboring subproblems
void TMOEAD::update_problem(CIndividual &indiv, int id)
{
	for (int i = 0; i<niche; i++)
	{
		int    k = population[id].table[i];
		double f1, f2;
		f1 = scalar_func(population[k].indiv.y_obj, population[k].namda,idealpoint, indivpoint);
		f2 = scalar_func(indiv.y_obj, population[k].namda,idealpoint, indivpoint);
		if (f2 < f1){
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
			//if (k < LSize)
			//	cout << "[" << mpi_rank << "] l update" << endl;
			//if (k >= pops_sub[mpi_rank] - RSize)
			//	cout << "[" << mpi_rank << "] r update" << endl;
		}
	}
}

bool TMOEAD::update_idealpoint_from_outside()
{
	bool isUpdateL = false;
	bool isUpdateR = false;
	if(numprocs >1){
		//只有多进程才使用，保证单进程可以执行
		int flag =0;
		MPI_Test(&req_idealpoint[L_RECV],&flag,&status_idealpoint);
		while(flag !=0)
		//if(flag != 0)
		{//接收到来自其他进程的理想点更新

			isUpdateL = update_idealpoint(recv_idealpoint[LEFT],LEFT);
			//cout<<"["<<mpi_rank<<"]("<<recv_idealpoint[0]<<","<<recv_idealpoint[1]<<") old("<<ij[0]<<","<<ij[1]<<")now("<<idealpoint[0]<<","<<idealpoint[1]<<")"<<endl;

			//继续等待接收更新
			MPI_Irecv(recv_idealpoint[LEFT],nobj,MPI_DOUBLE,prevRank,TAG_IDEALPOINT,MPI_COMM_WORLD,&req_idealpoint[L_RECV]);
			flag = 0;
			MPI_Test(&req_idealpoint[L_RECV],&flag,&status_idealpoint);
		}
		if (isUpdateL)
		{
			isIdealpointUpdate[LEFT] = true;
			isIdealpointUpdate[RIGHT] = false;
		}
		/*
		flag = 0;
		MPI_Test(&req_idealpoint[R_RECV], &flag, &status_idealpoint);
		while (flag != 0)
		{
			isUpdateR = update_idealpoint(recv_idealpoint[RIGHT],RIGHT);

			MPI_Irecv(recv_idealpoint[RIGHT],nobj,MPI_DOUBLE,nextRank,TAG_IDEALPOINT,MPI_COMM_WORLD,&req_idealpoint[R_RECV]);
			flag = 0;
			MPI_Test(&req_idealpoint[R_RECV],&flag,&status_idealpoint);
		}
		if (isUpdateR)
		{
			isIdealpointUpdate[RIGHT] = true;
			isIdealpointUpdate[LEFT] = false;
		}
		*/

	}
	return isUpdateL;//(isUpdateL || isUpdateR);
}

void TMOEAD::sendIdealpoint(double *buf,bool &isSend,int target,MPI_Request* request)
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
bool TMOEAD::update_idealpoint(CIndividual &ind)
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
	//isUpdate |= update_idealpoint_from_outside();
	if (isUpdate == true)
	{
		isIdealpointUpdate[LEFT] = true;
		isIdealpointUpdate[RIGHT] = true;
	}
	return isUpdate;
}


// update the reference point
bool TMOEAD::update_idealpoint(const double* y_obj,int side)
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
		isIdealpointUpdate[ side ] = true;
		isIdealpointUpdate[ 1 - side ] = true;
	}
	return isUpdate;
}

// recombination, mutation, update in MOEA/D
void TMOEAD::evolution()
{
	CIndividual child, child2;
	for(int i=LLSize; i<population.size()-RRSize; i++)
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

void  TMOEAD::execute(int mg, int irun, vector<double>& igd,vector<double>&hv, double &runtime,PTYPE ptype)
{
	//cout << "[" << mpi_rank << "]  ptype  " << ptype <<"   comm_i : "<< comm_i << endl;
	// mg: maximal number of generations
	init_uniformweight();
	init_neighborhood();

	if(numprocs >1){
		//只有多进程才使用，保证单进程也可以执行
		//通信邻居进程号
		prevRank = mpi_rank - 1;
		nextRank = mpi_rank + 1;
		if(mpi_rank == 0)
			prevRank = numprocs - 1;
		if(mpi_rank == (numprocs - 1))
			nextRank = 0;

		//分配缓冲空间
		//更新理想点使用的缓冲区
		if ((ptype == IDEAL_INDIV) || (ptype == IDEAL))
		{
			recv_idealpoint[LEFT] = new double[nobj];
			recv_idealpoint[RIGHT] = new double[nobj];
			send_idealpoint[LEFT] = new double[nobj];
			send_idealpoint[RIGHT] = new double[nobj];

			//非阻塞接收理想点更新
			//cout<<"["<<mpi_rank<<"] 1 before"<<endl;
			MPI_Irecv(recv_idealpoint[LEFT], nobj, MPI_DOUBLE, prevRank, TAG_IDEALPOINT, MPI_COMM_WORLD, &req_idealpoint[L_RECV]);
			MPI_Irecv(recv_idealpoint[RIGHT], nobj, MPI_DOUBLE, nextRank, TAG_IDEALPOINT, MPI_COMM_WORLD, &req_idealpoint[R_RECV]);

			//cout<<"["<<mpi_rank<<"] 1 after"<<endl;
			isBeginSend[LEFT] = false;
			isBeginSend[RIGHT] = false;
		}

		//个体更新接收
		//决策向量+目标向量+权重向量
		s_indiv =  nvar + nobj + 1;
		if ((ptype == IDEAL_INDIV) || (ptype == INDIV))
		{
			isLSend = false;
			isRSend = false;
			//cout<<"["<<mpi_rank<<"]LSize : "<<LSize<<endl;
			if (LLSize) {
				lrecv_indiv = new double[s_indiv * LLSize + 1];
				lsend_indiv = new double[s_indiv * LLSize + 1];
				//非阻塞等待接收左边邻居进程发送的更新个体
				//cout<<"["<<mpi_rank<<"] 2 before"<<endl;
				MPI_Irecv(lrecv_indiv, s_indiv*LLSize + 1, MPI_DOUBLE, prevRank, TAG_INDIV, MPI_COMM_WORLD, &req_indiv[L_RECV]);
				//cout<<"["<<mpi_rank<<"] 2 after"<<endl;
			}
			//cout<<"["<<mpi_rank<<"]RSize : "<<RSize<<endl;
			if (RRSize) {
				rrecv_indiv = new double[s_indiv*RRSize + 1];
				rsend_indiv = new double[s_indiv*RRSize + 1];
				//非阻塞等待接收右边邻居进程发送的更新个体
				//cout<<"["<<mpi_rank<<"] 3 before"<<endl;
				MPI_Irecv(rrecv_indiv, s_indiv * RRSize + 1, MPI_DOUBLE, nextRank, TAG_INDIV, MPI_COMM_WORLD, &req_indiv[R_RECV]);
				//cout<<"["<<mpi_rank<<"] 3 after"<<endl;
			}
		}
		//cout<<"["<<mpi_rank<<"] "<<endl;
		MPI_Barrier(MPI_COMM_WORLD);
	}

	//cout<<"["<<mpi_rank<<"] 4 before"<<endl;
	init_population();
	//cout<<"["<<mpi_rank<<"] 4 after"<<endl;

	if ((ptype == IDEAL_INDIV )|| (ptype == IDEAL))
	{
		if (numprocs > 1) {
			for (int i = 0; i < nobj; i++)
			{
				send_idealpoint[0][i] = idealpoint[i];
			}
			double *recv_utopian = new double[nobj*numprocs];

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

			MPI_Gatherv(send_idealpoint[0], nobj, MPI_DOUBLE, recv_utopian, rcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

			//统一理想点
			if (mpi_rank == 0)
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
			}
			//分发、同步理想点
			for (int i = 0; i < nobj; i++)
			{
				send_idealpoint[0][i] = idealpoint[i];
			}
			MPI_Bcast(send_idealpoint[0], nobj, MPI_DOUBLE, 0, MPI_COMM_WORLD);

			//更新进程内理想点
			for (int i = 0; i < nobj; i++) {
				idealpoint[i] = send_idealpoint[0][i];
			}
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
	if(IS_DEBUG){
		begintime = clock();
		//收集计算解集质量的必要数据,目标向量
		//等待，同步
		//MPI_Barrier(MPI_COMM_WORLD);
		if (mpi_rank == 0)
		{
			//载入真实前沿
			sprintf(filename, "%s/TestData/PF_Real/%s(%d).dat", WORKINGDIR, strTestInstance, nobj);
			loadpfront(filename, ps);
		}
		if(MODE_DEBUG)
		{
			gather_pop_y();
			//Root进程处理数据
			if(mpi_rank == 0)
			{
				//计算解集质量（IGD，HV)
				calc_qualities(igd_gen,hv_gen);
				igd.push_back(gen);  igd.push_back(igd_gen);
				hv.push_back(gen);  hv.push_back(hv_gen);
				cout<<"gen =  "<<gen<<"  hypervolume = "<<hv_gen<<"  "<<"  IGD = "<<igd_gen<<"  "<<endl;
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//等待，同步
		//MPI_Barrier(MPI_COMM_WORLD);
		endtime = clock();
		unevolvetime += endtime - begintime;
	}

	//开始执行进化
	for( gen=2; gen<=mg; gen++)
	{
		//cout<<"["<<mpi_rank<<"] gen"<<gen<<endl;
		evolution();

		if ((ptype == IDEAL_INDIV) || (ptype == IDEAL))
		{
			if (numprocs > 1)
			{
				if (isIdealpointUpdate[RIGHT])
					sendIdealpoint(send_idealpoint[LEFT], isBeginSend[LEFT], prevRank, &req_idealpoint[L_SEND]);
				if (isIdealpointUpdate[LEFT])
					sendIdealpoint(send_idealpoint[RIGHT], isBeginSend[RIGHT], nextRank, &req_idealpoint[R_SEND]);

				isIdealpointUpdate[LEFT] = false;
				isIdealpointUpdate[RIGHT] = false;
			}
			update_idealpoint_from_outside();
		}


		if ((ptype == IDEAL_INDIV) || (ptype == INDIV))
		{
			if (numprocs > 1 &&  (gen % comm_i == 0))
			{
				//更新个体
				if (LLSize)
				{
					//cout <<gen<< "  [" << mpi_rank << "] L unpdate inidv" << endl;
					if (isLSend)
					{
						MPI_Wait(&req_indiv[L_SEND], &status_indiv);
					}
					int m = 1;
					int lcount = 0;
					for (int p = 0; p < LLSize; p++)
					{
						if (recordUpdateofIndiv[p])
						{

							TSOP &targetP = population[p];

							lsend_indiv[m++] = targetP.array[0];
							for (int i = 0; i < nvar; i++)
								lsend_indiv[m++] = targetP.indiv.x_var[i];

							for (int i = 0; i < nobj; i++)
								lsend_indiv[m++] = targetP.indiv.y_obj[i];

							//for (int i = 0; i < nobj; i++)
							//	lsend_indiv[m++] = targetP.namda[i];


							lcount++;
						}
					}
					if (lcount > 0)
					{
						lsend_indiv[0] = (double)lcount;
						//cout << gen << "  [" << mpi_rank << "] -- > " << prevRank << "  update indiv : " << lcount << endl;
						//cout<<"["<<mpi_rank<<"]"<<size_node<<"\t"<<s_indiv-1-1-2<<endl;
						MPI_Isend(lsend_indiv, LLSize * s_indiv + 1, MPI_DOUBLE, prevRank, TAG_INDIV, MPI_COMM_WORLD, &req_indiv[L_SEND]);
						//发送指令完毕，非阻塞，继续执行接下来计算,边传输边计算
						//cout<<"["<<mpi_rank<<"] send"<<endl;

						isLSend = true;
					}
				}
				//cout<<"["<<mpi_rank<<"] 4"<<endl;
				if (RRSize)
				{
					//cout <<gen <<  "  [" << mpi_rank << "] R unpdate inidv" << endl;
					if (isRSend)
					{
						MPI_Wait(&req_indiv[R_SEND], &status_indiv);
					}
					int m = 1;
					int rcount = 0;
					for (int p = pops_sub[mpi_rank] - RRSize; p < pops_sub[mpi_rank]; p++)
					{
						if (recordUpdateofIndiv[p])
						{
							TSOP &targetP = population[p];

							rsend_indiv[m++] = targetP.array[0];
							for (int i = 0; i < nvar; i++)
								rsend_indiv[m++] = targetP.indiv.x_var[i];
							for (int i = 0; i < nobj; i++)
								rsend_indiv[m++] = targetP.indiv.y_obj[i];
							//for (int i = 0; i < nobj; i++)
							//	rsend_indiv[m++] = targetP.namda[i];



							rcount++;
						}
					}
					if (rcount > 0)
					{
						rsend_indiv[0] = (double)rcount;
						//cout << gen << "  [" << mpi_rank << "] -- > " << nextRank << "  update indiv : " << rcount << endl;
						//cout<<"["<<mpi_rank<<"]"<<size_node<<"\t"<<s_indiv-1-1-2<<endl;
						MPI_Isend(rsend_indiv, RRSize * s_indiv + 1, MPI_DOUBLE, nextRank, TAG_INDIV, MPI_COMM_WORLD, &req_indiv[R_SEND]);
						//发送指令完毕，非阻塞，继续执行接下来计算,边传输边计算
						//cout<<"["<<mpi_rank<<"] send"<<endl;

						isRSend = true;
					}
				}
				//cout<<"["<<mpi_rank<<"] 5"<<endl;

				//重置更新标志
				for (int i = 0; i < pops_sub[mpi_rank]; i++)
					recordUpdateofIndiv[i] = false;
			}

			update_indiv_from_outside();
		}


		if(IS_DEBUG && (gen % I_OUT == 0)){


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
					if (mpi_rank == 0)
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
		if ((ptype == IDEAL_INDIV) || (ptype == INDIV))
		{
			//清理操作，完成操作
			isDone[L_RECV] = isDone[R_RECV] = 1;
			if (LLSize) {
				isDone[L_RECV] = 0;
				MPI_Irecv(&b_isdone[LEFT], 1, MPI_INT, prevRank, TAG_DONE, MPI_COMM_WORLD, &req_done[L_RECV]);
			}
			if (RRSize) {
				isDone[R_RECV] = 0;
				MPI_Irecv(&b_isdone[RIGHT], 1, MPI_INT, nextRank, TAG_DONE, MPI_COMM_WORLD, &req_done[R_RECV]);
			}

			int flag = 1;
			if (isLSend)
			{
				MPI_Test(&req_indiv[L_SEND], &flag, &status_indiv);
				if (flag)
					isLSend = false;
			}
			if (isLSend)
				isDone[L_SEND] = 0;
			else
			{
				isDone[L_SEND] = 1;
				if (LLSize)
					MPI_Isend(&isDone[L_SEND], 1, MPI_INT, prevRank, TAG_DONE, MPI_COMM_WORLD, &req_done[L_SEND]);
			}


			if (isRSend)
			{
				MPI_Test(&req_indiv[R_SEND], &flag, &status_indiv);
				if (flag)
					isRSend = false;
			}

			if (isRSend)
				isDone[R_SEND] = 0;
			else
			{
				isDone[R_SEND] = 1;
				if (RRSize)
					MPI_Isend(&isDone[R_SEND], 1, MPI_INT, nextRank, TAG_DONE, MPI_COMM_WORLD, &req_done[R_SEND]);
			}

			while (1)
			{
				if (!isDone[L_SEND])
				{
					if (isLSend)
					{
						MPI_Test(&req_indiv[L_SEND], &flag, &status_indiv);
						if (flag)
							isLSend = false;
					}
					if (!isLSend)
					{
						isDone[L_SEND] = 1;
						MPI_Isend(&isDone[L_SEND], 1, MPI_INT, prevRank, TAG_DONE, MPI_COMM_WORLD, &req_done[L_SEND]);
					}
				}
				if (!isDone[R_SEND])
				{
					if (isRSend)
					{
						MPI_Test(&req_indiv[R_SEND], &flag, &status_indiv);
						if (flag)
							isRSend = false;
					}
					if (!isRSend)
					{
						isDone[R_SEND] = 1;
						MPI_Isend(&isDone[R_SEND], 1, MPI_INT, nextRank, TAG_DONE, MPI_COMM_WORLD, &req_done[R_SEND]);
					}
				}
				if (!isDone[L_RECV]) {
					MPI_Test(&req_done[L_RECV], &flag, &status_done);
					if (flag && b_isdone[LEFT])
					{//左边完成进化
						isDone[L_RECV] = 1;
						//	MPI_Cancel(&req_indiv[L_RECV]);
						//MPI_Request_free(&req_indiv[L_RECV]);
					}
					else
					{
						if (ptype == IDEAL_INDIV)
							update_idealpoint_from_outside();
						if (LLSize)
							getAndUpdateIndivs(lrecv_indiv, s_indiv*LLSize + 1, prevRank, &req_indiv[L_RECV]);

					}
				}
				if (!isDone[R_RECV])
				{
					MPI_Test(&req_done[R_RECV], &flag, &status_done);
					if (flag && b_isdone[RIGHT])
					{//右边完成进化
						isDone[R_RECV] = 1;
						//	MPI_Cancel(&req_indiv[R_RECV]);
						//	MPI_Request_free(&req_indiv[R_RECV]);
					}
					else
					{
						if (ptype == IDEAL_INDIV)
							update_idealpoint_from_outside();
						if (RRSize)
							getAndUpdateIndivs(rrecv_indiv, s_indiv*RRSize + 1, nextRank, &req_indiv[R_RECV]);
					}
				}
				if (isDone[L_SEND] && isDone[L_RECV] && isDone[R_SEND] && isDone[R_RECV])
					break;
			}
		}

		MPI_Barrier(MPI_COMM_WORLD);
		if ((ptype == IDEAL_INDIV) || (ptype == IDEAL))
		{
			update_idealpoint_from_outside();
			MPI_Cancel(&req_idealpoint[L_RECV]);
			MPI_Request_free(&req_idealpoint[L_RECV]);
			MPI_Cancel(&req_idealpoint[R_RECV]);
			MPI_Request_free(&req_idealpoint[R_RECV]);


			//等待，同步
			if (isBeginSend[LEFT])
			{
				int flag = 1;
				MPI_Test(&req_idealpoint[L_SEND], &flag, &status_idealpoint);
				if (flag == 0)
				{
					MPI_Cancel(&req_idealpoint[L_SEND]);
					MPI_Request_free(&req_idealpoint[L_SEND]);
				}
			}
			if (isBeginSend[RIGHT])
			{
				int flag = 1;
				MPI_Test(&req_idealpoint[R_SEND], &flag, &status_idealpoint);
				if (flag == 0)
				{
					MPI_Cancel(&req_idealpoint[R_SEND]);
					MPI_Request_free(&req_idealpoint[R_SEND]);
				}
			}
		}

		if ((ptype == IDEAL_INDIV) || (ptype == INDIV))
		{
			if (LLSize)
			{
				MPI_Cancel(&req_indiv[L_RECV]);
				MPI_Request_free(&req_indiv[L_RECV]);
			}

			if (isLSend)
			{
				int flag = 1;
				MPI_Test(&req_indiv[L_SEND], &flag, &status_indiv);
				if (flag == 0)
				{
					MPI_Cancel(&req_indiv[L_SEND]);
					MPI_Request_free(&req_indiv[L_SEND]);
				}
			}
			if (RRSize)
			{
				MPI_Cancel(&req_indiv[R_RECV]);
				MPI_Request_free(&req_indiv[R_RECV]);
			}
			if (isRSend)
			{
				int flag = 1;
				MPI_Test(&req_indiv[R_SEND], &flag, &status_indiv);
				if (flag == 0)
				{
					MPI_Cancel(&req_indiv[R_SEND]);
					MPI_Request_free(&req_indiv[R_SEND]);
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
		sprintf(filename, "%s/TestData/POP/POP_%s_%s(%d)_%d_%dR%d",WORKINGDIR,STR_MOEAD, strTestInstance,nobj, popsize, total_gen, irun);
		sprintf(filename, "%s_%dnp_%s_%s",filename,numprocs,getElemTypeName(paraType),getElemTypeName(paraTopology) );
		sprintf(filename, "%s.dat",filename);
		save_population_front_split(pop_y,pop_index,index_split,filename);//保存扇形存档种群的非劣前沿
		pop_index.clear();
		//过滤取出前沿
		vector <vector<double> >  pf_y;
		vector<int> pf_index;
		population2front(pop_y,pf_y,pf_index);
		//保存整体前沿
		sprintf(filename, "%s/TestData/POF/PF_%s_%s(%d)_%d_%dR%d",WORKINGDIR,STR_MOEAD, strTestInstance,nobj, popsize, total_gen, irun);
		sprintf(filename, "%s_%dnp_%s_%s",filename,numprocs,getElemTypeName(paraType),getElemTypeName(paraTopology) );
		sprintf(filename, "%s.dat",filename);
		save_population_front_split(pf_y,pf_index,index_split,filename);//保存扇形存档种群的非劣前沿
		//保存决策变量
		sprintf(filename, "%s/TestData/POS/POS_%s_%s(%d)_%d_%dR%d",WORKINGDIR,STR_MOEAD, strTestInstance,nobj, popsize, total_gen, irun);
		sprintf(filename, "%s_%dnp_%s_%s",filename,numprocs,getElemTypeName(paraType),getElemTypeName(paraTopology) );
		sprintf(filename, "%s.dat",filename);
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
		if ((ptype == IDEAL_INDIV) || (ptype == IDEAL))
		{
			delete[]recv_idealpoint[LEFT];
			recv_idealpoint[LEFT] = NULL;
			delete[]recv_idealpoint[RIGHT];
			recv_idealpoint[RIGHT] = NULL;
			delete[]send_idealpoint[LEFT];
			send_idealpoint[LEFT] = NULL;
			delete[]send_idealpoint[RIGHT];
			send_idealpoint[RIGHT] = NULL;
		}
		/*for(int i=0;i<nobj;i++)
		{
		delete []recv_anchorpoint[LEFT][i];
		recv_anchorpoint[LEFT][i] = NULL;
		delete []recv_anchorpoint[RIGHT][i];
		recv_anchorpoint[RIGHT][i] = NULL;
		delete []send_anchorpoint[LEFT][i];
		send_anchorpoint[LEFT][i] = NULL;
		delete []send_anchorpoint[RIGHT][i];
		send_anchorpoint[RIGHT][i] = NULL;
		}*/

		if ((ptype == IDEAL_INDIV) || (ptype == INDIV))
		{
			if (LLSize != 0)
			{
				delete[]lrecv_indiv;
				lrecv_indiv = NULL;
				delete[]lsend_indiv;
				lsend_indiv = NULL;
			}
			if (RRSize != 0)
			{
				delete[]rrecv_indiv;
				rrecv_indiv = NULL;
				delete[]rsend_indiv;
				rsend_indiv = NULL;
			}
		}
		//cout << "[" << mpi_rank << "]finish" << endl;
	}
}


void TMOEAD::save_front(char saveFilename[1024])
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

void TMOEAD::save_pos(char saveFilename[1024])
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


void TMOEAD::operator=(const TMOEAD &emo)
{
    popsize        = emo.popsize;
	population  = emo.population;
	indivpoint  = emo.indivpoint;
	niche       = emo.niche;
}

void TMOEAD::gather_pop_y()
{
	//发送各个个体的目标向量，用于计算IGD和HV指标
	int m=0;
	for(int i =LLSize;i<LLSize + pops_sub_self[mpi_rank];i++)
	{
	 //for(int j =0;j<nvar;j++)
		//b_send[m++] = population[i].indiv.x_var[j];
		for(int j=0;j<nobj;j++){
			b_send[m++] = population[i].indiv.y_obj[j];
		}
	}
	int size_send = nobj * pops_sub_self[mpi_rank];

	//收集数据
	int *rcounts = new int[numprocs];
	for(int i=0;i<numprocs;i++)
		rcounts[i] = nobj*pops_sub_self[i];
	int *displs = new int[numprocs];
	displs[0] = 0;
	for(int i=1;i<numprocs;i++)
		displs[i] = displs[i-1] + rcounts[i-1];

	MPI_Gatherv(b_send,size_send,MPI_DOUBLE,b_recv,rcounts,displs,MPI_DOUBLE,0,MPI_COMM_WORLD);
	delete rcounts;
	rcounts = NULL;
	delete displs;
	displs = NULL;
}

void TMOEAD::gather_populations()
{
	//cout << "["<<mpi_rank<<"]mark" << endl;
	//发送各个个体的目标向量，用于计算IGD和HV指标
	int m = 0;
	for (int i = LLSize; i < LLSize + pops_sub_self[mpi_rank]; i++)
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
	MPI_Gatherv(b_send,s_send,MPI_DOUBLE,b_recv,rcounts,displs,MPI_DOUBLE,0,MPI_COMM_WORLD);
	delete rcounts;
	rcounts = NULL;
	delete displs;
	displs = NULL;
	//cout << "["<<mpi_rank<<"]mark end" << endl;
}

void TMOEAD::calc_qualities(double &igd,double& hv)
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

double TMOEAD::calc_distance(vector<vector<double> >& ParetoFront)
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

 double TMOEAD::calc_hypervolume(vector<vector<double> >& ParetoFront)
{
	wfg hv;
	double result = hv.compute(ParetoFront,hv_ref_ub);

	return result / cubeHV;
}

void TMOEAD::getAndUpdateIndivs(double* b_recv_indiv,int datasize,int source,MPI_Request *request)
{
	//if(mpi_rank == 0)
		//cout<<"["<<mpi_rank<<"] LSize : "<<LSize<<"\tRSize : "<<RSize<<endl;
	int flag =0;
	MPI_Test(request,&flag,&status_indiv);
	if(flag !=0)
	//if(flag !=0)
	{//接收到来自进程的个体更新
		//从缓存中获取数据
		int m = 0;
		int count = (int)b_recv_indiv[m++];
		//cout << "[" << mpi_rank << "] receive indiv from  " << source << "   size:  " << count << endl;
		for (int p = 0; p < count; p++)
		{
			int index = floor(b_recv_indiv[m++] + 0.5) - pops_sub_start[mpi_rank];
			CIndividual indiv;
			for (int i = 0; i < nvar; i++)
				indiv.x_var[i] = b_recv_indiv[m++];
			for (int i = 0; i < nobj; i++)
				indiv.y_obj[i] = b_recv_indiv[m++];


			//update idealpoint
			update_idealpoint(indiv);

			//根据权重向量计算对应种群下标
			//int index = floor(lambda[0] * n_H + 0.5) - pops_sub_start[mpi_rank];
				//cout<<"["<<mpi_rank<<"]"<<index<<endl;
			for (int i = 0; i < niche; i++)
			{
				int k = population[index].table[i];

				double f1, f2;
				f1 = scalar_func(population[k].indiv.y_obj, population[k].namda,idealpoint, indivpoint);
				f2 = scalar_func(indiv.y_obj, population[k].namda, idealpoint,indivpoint);
				if (f2 < f1) {
					population[k].indiv = indiv;
					recordUpdateofIndiv[k] = false;
				}
			}
			//lambda.clear();
		}

		MPI_Irecv(b_recv_indiv, datasize, MPI_DOUBLE, source, TAG_INDIV, MPI_COMM_WORLD, request);
	}
}

void TMOEAD::sendIndiv(double* b_send_indiv,int index,int target,MPI_Request *request,bool &isBegin)
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

void TMOEAD::update_indiv_from_outside()
{
	if (numprocs > 1)
	{
		//左边更新
		if (LLSize)
			getAndUpdateIndivs(lrecv_indiv, s_indiv*LLSize + 1, prevRank, &req_indiv[L_RECV]);

		//右边更新
		if (RRSize)
			getAndUpdateIndivs(rrecv_indiv, s_indiv*RRSize + 1, nextRank, &req_indiv[R_RECV]);
	}

}

void population2front(vector <TSOP>  &mypopulation, vector <CIndividual>  &population_front)
{
	vector<int> nDominated;
	for (int n = 0; n<mypopulation.size(); n++)
		nDominated.push_back(0);


	for (int k = 0; k<mypopulation.size(); k++)
		for (int j = k + 1; j<mypopulation.size(); j++)
		{
			TCompare tresult = ParetoCompare(mypopulation[k].indiv.y_obj, mypopulation[j].indiv.y_obj);
			if (tresult == _Pareto_Dominated)
				nDominated[k]++;
			else if (tresult == _Pareto_Dominating)
				nDominated[j]++;
		}

	for (int n = 0; n<mypopulation.size(); n++)
		if (nDominated[n] == 0)
			population_front.push_back(mypopulation[n].indiv);

	nDominated.clear();
}
void population2front(vector <TSOP>  &mypopulation, vector <vector<double> >  &population_front)
{
	vector<int> nDominated;
	for (int n = 0; n<mypopulation.size(); n++)
		nDominated.push_back(0);


	for (int k = 0; k<mypopulation.size(); k++)
		for (int j = k + 1; j<mypopulation.size(); j++)
		{
			TCompare tresult = ParetoCompare(mypopulation[k].indiv.y_obj, mypopulation[j].indiv.y_obj);
			if (tresult == _Pareto_Dominated)
				nDominated[k]++;
			else if (tresult == _Pareto_Dominating)
				nDominated[j]++;
		}

	for (int n = 0; n<mypopulation.size(); n++)
		if (nDominated[n] == 0)
			population_front.push_back(mypopulation[n].indiv.y_obj);

	nDominated.clear();
}
void save_population(vector <CIndividual>  &mypopulation, char saveFilename[1024])
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
