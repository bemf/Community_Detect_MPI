#ifndef __CEC09_H_
//#include "cec09.h"
#endif

#include "../global.h"
#include "../recomb.h"
#include "../common.h"

#include "cheaind.h"

#include "cheaalg.h"
#include "../TestFunction/testInstance.h"
#include<iostream>
using namespace std;
CHEA::CHEA(int pop_size, int hype_intercept){

	start = clock();
	popsize = pop_size;
	hyperplane_Intercept = hype_intercept;
	population.reserve(popsize);
	TrueNadirPoint.reserve(nobj);
	ReferencePoint.reserve(nobj);
	IdealPoint.reserve(nobj);
	per_intercept = 1.0 / hyperplane_Intercept;

	size_subp_on_edge = 0;
	maxLayer = 1;

	//����ֲ�������������С
	pops_sub_self = new int[numprocs];
	int mod_tmp = popsize %numprocs;
	int v_tmp = popsize / numprocs;
	if (mod_tmp == 0)
		pops_sub_self[0] = v_tmp;
	else
		pops_sub_self[0] = v_tmp + 1;

	for (int i = 1; i<numprocs; i++)
		pops_sub_self[i] = v_tmp + (mod_tmp > i ? 1 : 0);

	//cout << "[" << mpi_rank << "] " << pops_sub_self[mpi_rank] <<endl;

	//����Ǩ�ƾ�Ӣ���建��ռ�
	int tmp = nvar + nobj;
    s_send = tmp * pops_sub_self[mpi_rank];
    s_recv = tmp * popsize;
    b_send = new double[s_send];

    cs_send = tmp * popsize;
    cs_recv = tmp * popsize*numprocs;
	c_send=new double[cs_send];

	if (mpi_rank == 0)
    {
		b_recv = new double[s_recv];
        c_recv=new double[cs_recv];
    }
	// initialize ideal point
	for (int n = 0; n < nobj; n++)
		IdealPoint.push_back(1.0e+30);

    //send_indivs.resize(numprocs);
}

CHEA::~CHEA(){
	delete []b_send;
	if (mpi_rank == 0)
		delete []b_recv;
	lambda.clear();
	population.clear();
	ps.clear();
	TrueNadirPoint.clear();
	IdealPoint.clear();
	ReferencePoint.clear();

    delete []c_send;
    if(mpi_rank==0)delete[]c_recv;
}


//��ʼ��
void CHEA::init_population()
{

	population[0].indiv.rnd_init();
	//��ʼ����ֵ

	for (int j = 0; j<nobj; j++)
	{
		TrueNadirPoint.push_back(population[0].indiv.y_obj[j]);
		IdealPoint.push_back(population[0].indiv.y_obj[j]);
		ReferencePoint.push_back(IdealPoint[j] + 1e3 * (TrueNadirPoint[j] - IdealPoint[j]));
	}
	for (int n = 1; n < popsize; n++)
	{
		population[n].indiv.rnd_init();
		update_extreme_point(population[n].indiv);
	}

	//�����ʼ�����Ӧ������ı��
	update_partition();
}

//�ֽ������⡢������۲�����
void CHEA::gen_uniformweight(int start_obj_index, int max_value_left, vector<int> coordinate, int H)
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
		gen_uniformweight(start_obj_index - 1, max_value_left - i, coordinate, H);
	}
}

// initialize a set of evely-distributed weight vectors
void CHEA::init_uniformweight()
{
	vector<int> coordinate(nobj, 0);
	gen_uniformweight(nobj - 1, n_H, coordinate, n_H);

	//������Ȩ��
	for (int i = nobj-1; i >= 0; i--)
	{
		coordinate[i] = hyperplane_Intercept;
		subp_index_on_edge.push_back(cal_index_from_lambda(coordinate, hyperplane_Intercept));
		coordinate[i] = 0;
	}
	size_subp_on_edge = subp_index_on_edge.size();
	coordinate.clear();

	for (int i = 0; i < popsize; i++)
	{
		CHEA_SOP sop;
		sop.sectorialindex = i;
		sop.V_obj_o = lambda[i];
		population.push_back(sop);
	}

	//Ϊÿ�����̷���Ȩ������
	allocateForProcess();
	//allocateForProcess_Linear();

	if (numprocs > 1)
	{
		//�������˽ṹ
		//����
		if (paraTopology == UNI_RING)
		{
			neighbor_partitions.push_back((mpi_rank + 1 + numprocs) % numprocs);
			partition_niche = neighbor_partitions.size();
		}
		else if (paraTopology == BI_RING)
		{
			//˫��
			neighbor_partitions.push_back((mpi_rank + 1 + numprocs) % numprocs);
			neighbor_partitions.push_back((mpi_rank - 1 + numprocs) % numprocs);
			partition_niche = neighbor_partitions.size();
		}
		else if (paraTopology == FULL)
		{
			//��ȫ����
			for (int i = 0; i < numprocs; i++)
			{
				if (mpi_rank != i)
					neighbor_partitions.push_back(i);
			}
			partition_niche = neighbor_partitions.size();
		}
		else if (paraTopology == NEI)
		{
			//�����ھӹ�ϵ������
			//���串��������Ȩ��
			partition_overlap(maxLayer);
			//cout << "[" << mpi_rank << "] w2" << endl;
			//������������ھ�
			partition_neighbor();
		}

	}
}



bool CHEA::lambda_compare(vector<int> &l1, vector<int> &l2)
{
	int dim = l1.size();
	for (int i = 0; i < dim - 1; i++)
	{
		if (l1[i] < l2[i])
			return true;
		else if (l1[i] > l2[i])
			return false;
	}
	return false;
}

void CHEA::allocateForProcess()
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
						int index = j++;
						int l_dim = lambda[index][nobj - 1];
						while (j < popsize && lambda[j][nobj - 1] == l_dim)
						{
							if (isAllocated[j] == false && lambda_compare(lambda[j], lambda[index]))
								index = j;
							j++;
						}
						inQueue.push(index);
						inSet.insert(index);
						break;
					}
			}
			int index = inQueue.front();
			inQueue.pop();
			//����ӳ��
			if (isAllocated[index] == false)
			{
				//Ȩ�ع���partition
				sop_partitions[index].push_back(i);
				partition_indexs[i].insert(index);

				isAllocated[index] = true;
				count++;
				if (low > lambda[index][nobj - 1])
					low = lambda[index][nobj - 1];
				//�����ھӼ���ѡ�����
				vector<int> vicinities = gen_vicinities(lambda[index],hyperplane_Intercept);
				//cout << "size : " << vicinities.size() << endl;
				for (int j = 0; j < vicinities.size(); j++)
					if (isAllocated[vicinities[j]] == false && inSet.count(vicinities[j]) == 0)
					{
						inQueue.push(vicinities[j]);
						inSet.insert(vicinities[j]);
					}
			}
		}
		//һ�����̷������
		inSet.clear();
		if ((low > 0) && !inQueue.empty())
		{
			int index = inQueue.front();
			inQueue.pop();
			while (!inQueue.empty())
			{
				int index2 = inQueue.front();
				inQueue.pop();

				/*
				if (lambda[index2][0] < lambda[index][0])
				index = index2;
				else if ((lambda[index2][0] == lambda[index][0]) && (lambda[index2][nobj - 1] > lambda[index][nobj - 1]))
				index = index2;
				*/
				if (lambda_compare(lambda[index2], lambda[index]))
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


void CHEA::allocateForProcess_Linear()
{
	sop_partitions = vector<vector<int> >(popsize);
	partition_indexs = vector<set<int> >(numprocs);

	int index = 0;
	for (int i = 0; i < numprocs; i++)
	{
		int count = 0;
		while (count < pops_sub_self[i])
		{
			//Ȩ�ع���partition
			sop_partitions[index].push_back(i);
			partition_indexs[i].insert(index++);
			count++;
		}
	}

}


void CHEA::partition_overlap(int maxLayer)
{
	partition_overlap_indexs = vector<set<int> >(numprocs);
	for (int i = 0; i < numprocs; i++)
	{
		set<int> lastLay(partition_indexs[i]);

		for (int l = 0; l < maxLayer; l++)
		{
			set<int> newLay;
			set<int>::iterator it = lastLay.begin();
			set<int>::iterator it_end = lastLay.end();
			while (it != it_end)
			{
				vector<int> vicinities = gen_vicinities(lambda[*it],hyperplane_Intercept);
				for (int k = 0; k < vicinities.size(); k++)
				{
					if (partition_indexs[i].count(vicinities[k]) == 0 && partition_overlap_indexs[i].count(vicinities[k]) == 0)
					{
						partition_overlap_indexs[i].insert(vicinities[k]);
						newLay.insert(vicinities[k]);
						sop_partitions[vicinities[k]].push_back(i);
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

void CHEA::partition_neighbor()
{
	//cout << "[" << mpi_rank << "] : " << sop_partitions.size() << endl;
	//vector<vector<int> > nicheC(numprocs,vector<int>(numprocs,0));
	for (int i = 0; i < sop_partitions.size(); i++)
	{
		//cout << "[" << mpi_rank << "] sop : " << sop_partitions[i].size() << endl;
		//sop_partitions[i].size() > 1//�й��ã����ǻ��า������
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

//���·�����嵽��Ӧ����������
void CHEA::update_partition()
{

	//cout << "[" << mpi_rank << "] update_partition " << endl;
	//������Ⱥ��ÿ������ı��ֵ
	int n;
	for (n = 0; n < popsize; n++)
    {
		population[n].indiv.obj_index(IdealPoint, hyperplane_Intercept);
    }
	//����ÿ��������������Ķ�Ӧ��ϵ

	vector <bool>  sectorpopFlag(popsize,false);

	vector <CHEAInd>  initial_rest_indi_pop;
	for (int n = 0; n< popsize; n++)//
	{
		if (population[n].sectorialindex == population[n].indiv.sectorialindex)
		{
			sectorpopFlag[population[n].sectorialindex] = true;
		}
		else{
			CHEAInd replacedInd = population[n].indiv;
			bool result;
			//while(true){
				sectorpopFlag[replacedInd.sectorialindex] = true;
				result = Pareto_HyperVolume_compare_sectorialgrid(replacedInd, replacedInd);
				//if (result == false)
				//	break;
			//}
			initial_rest_indi_pop.push_back(replacedInd);
		}
	}

	int size_of_rest_indi = initial_rest_indi_pop.size();
	initial_rest_indi_pop.reserve(size_of_rest_indi);
	vector<vector<int> > V_obj_rest_indi;
	V_obj_rest_indi.reserve(size_of_rest_indi);

	vector<int>cal_V_obj(nobj,0);
	for (int i = 0; i < size_of_rest_indi; i++)
	{
		initial_rest_indi_pop[i].cal_V_obj(IdealPoint, hyperplane_Intercept, cal_V_obj);
		V_obj_rest_indi.push_back(cal_V_obj);
	}
	//cout << "[" << mpi_rank << "] restore left subproblems" << endl;
	//���ʣ��ĸ���
	for (int i = 0; i < popsize; i++)
	{
		if (sectorpopFlag[i] == false)
		{
			CHEA_SOP& subproblem = population[i];

			int min_index_dist = 0;
			for (int j = 0; j < nobj; j++)
			{
				min_index_dist = pow(V_obj_rest_indi[0][j] - subproblem.V_obj_o[j], 2.0);
			}

			int min_diff_index = 0;
			int rest_size = initial_rest_indi_pop.size();

			for (int k = 1; k < rest_size; k++)
			{
				int index_dist = 0;
				for (int j = 0; j < nobj; j++)
				{
					index_dist = pow(V_obj_rest_indi[k][j] - subproblem.V_obj_o[j], 2.0);
				}

				if (index_dist < min_index_dist || (index_dist == min_index_dist && rnd_uni(&rnd_uni_init) > 0.5))
				{
					min_index_dist = index_dist;
					min_diff_index = k;
				}
			}


			population[i].indiv = initial_rest_indi_pop[min_diff_index];
			sectorpopFlag[i] = true;
			//rest_ind = initial_rest_indi_pop[min_diff_index];
			initial_rest_indi_pop[min_diff_index] = initial_rest_indi_pop[rest_size - 1];
			V_obj_rest_indi[min_diff_index] = V_obj_rest_indi[rest_size - 1];
			//initial_rest_indi_pop[rest_size - 1] = rest_ind;
			initial_rest_indi_pop.pop_back();
			V_obj_rest_indi.pop_back();
		}
	}
	initial_rest_indi_pop.clear();
	V_obj_rest_indi.clear();
	sectorpopFlag.clear();
	//cout << "[" << mpi_rank << "] update partition finish" << endl;

}



//����������ȫ�ֲο���
bool CHEA::update_ideal_point(const double *y_obj)
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
		if (bUpdate_item )
		{
			ReferencePoint[j] = (TrueNadirPoint[j]) + 1e3 *(TrueNadirPoint[j] - IdealPoint[j]);
			bUpdate_item = false;
		}
	}

	isIdealpointUpdate |= bIdealUpdated;
	return bIdealUpdated;
}

bool CHEA::update_extreme_point(CHEAInd& ind)
{
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

//ѡ����������ڱ�����ָ��Ľ�����ѡ������
int  CHEA::tour_selection_hv(vector <CHEA_SOP>  &population)
{
	//int p1 = int(rnd_uni(&rnd_uni_init)*population.size());
	//int p2 = int(rnd_uni(&rnd_uni_init)*population.size());
	int c1 = int(rnd_uni(&rnd_uni_init)*(pops_sub_self[mpi_rank] -1));
	int c2 = int(rnd_uni(&rnd_uni_init)*(pops_sub_self[mpi_rank] - c1));

	set<int>::iterator it = partition_indexs[mpi_rank].begin();

	for (int i = 0; i < c1; i++)
		it++;
	int p1 = *it;
	it++;
	for (int i = c1+1; i < c2; i++)
		it++;
	int p2 = *it;

	//cout << "[" << mpi_rank << "] sp1 = " << p1 << " ; cp2 = " << p2 << endl;

	double hv1 = tour_selection_hv_difference(p1, population);
	double hv2 = tour_selection_hv_difference(p2, population);

	//cout << "[" << mpi_rank << "]  hv1 = " << hv1 << "  ; hv2 = " << hv2 << endl;
	if (hv1 >= hv2)//ѡ������ֵ�ϴ�ĸ���
		return p1;
	else
		return p2;
}
int  CHEA::tour_selection_Q(vector <CHEA_SOP>  &population,MGraph *G)
{
	//int p1 = int(rnd_uni(&rnd_uni_init)*population.size());
	//int p2 = int(rnd_uni(&rnd_uni_init)*population.size());
	int c1 = int(rnd_uni(&rnd_uni_init)*(pops_sub_self[mpi_rank] -1));
	int c2 = int(rnd_uni(&rnd_uni_init)*(pops_sub_self[mpi_rank] - c1));

	set<int>::iterator it = partition_indexs[mpi_rank].begin();

	for (int i = 0; i < c1; i++)
		it++;
	int p1 = *it;
	it++;
	for (int i = c1+1; i < c2; i++)
		it++;
	int p2 = *it;

	//cout << "[" << mpi_rank << "] sp1 = " << p1 << " ; cp2 = " << p2 << endl;
    double qmetric1=0.0,qmetric2=0.0;
    //qmetric1=calc_q(p1,population,G);
    //qmetric2=calc_q(p2,population,G);

    qmetric1=-population[p1].indiv.y_obj[0]-population[p1].indiv.y_obj[1];
    qmetric2=-population[p2].indiv.y_obj[0]-population[p2].indiv.y_obj[1];
    //qmetric1=-population[p1].indiv.y_obj[0]*population[p1].V_obj_o[0]-population[p1].indiv.y_obj[1]*population[p1].V_obj_o[1];
    //qmetric2=-population[p2].indiv.y_obj[0]*population[p2].V_obj_o[0]-population[p2].indiv.y_obj[1]*population[p2].V_obj_o[1];



	//cout << "[" << mpi_rank << "]  hv1 = " << hv1 << "  ; hv2 = " << hv2 << endl;
	if (qmetric1 >= qmetric2)//ѡ������ֵ�ϴ�ĸ���
		return p1;
	else
		return p2;
}

double CHEA::calc_q(int index,vector <CHEA_SOP>  &population,MGraph *G)
{
     double qmetric=0.0;
    vector<vector<double> >model;
    getModel(population[index].indiv.x_var,model);

    int model_size=model.size();
    for(int j=0; j < model_size; j++ )
    {
     double q1=(getL_NRA(model[j],model[j],G)/2)/G->numEdges;
     int kDegree=getDegree(model[j],G);
     double q2=(float)(kDegree*kDegree)/(4*G->numEdges*G->numEdges);

     qmetric+=q1-q2;
    }
    return qmetric;
}


//�Ƚϸ���������Ľ�
bool CHEA::Pareto_HyperVolume_compare_sectorialgrid(const CHEAInd& ind ,CHEAInd& replacedInd)
{
	bool bReplaced = false;

	double contributing1, contributing2;
    //printf("sectorialindex:%d,obj:%f\n", ind.sectorialindex,ind.y_obj[0]);
	CHEA_SOP &subproblem = population[ind.sectorialindex];
	subproblem.indiv.obj_index(IdealPoint,hyperplane_Intercept);
	if (subproblem.sectorialindex == subproblem.indiv.sectorialindex)
	{
		//�ڴ˴����ö�̬�ο�����߾�̬�ο���

		vector<double> ref_cal(nobj,0);

		subproblem.indiv.cal_k_value(IdealPoint,hyperplane_Intercept);

		double k = ind.k_value > subproblem.indiv.k_value ? ind.k_value : subproblem.indiv.k_value;

		for (int i = 0; i < nobj; i++)
		{
			ref_cal[i] = (IdealPoint[i] + k * per_intercept * (subproblem.V_obj_o[i] + 1));
		}

		contributing1 = GetHyperVolume(ind, ref_cal);//ReferencePoint);//
		contributing2 = GetHyperVolume(subproblem.indiv, ref_cal);//ReferencePoint);//

		ref_cal.clear();

		if (contributing1 > contributing2)
		{
			replacedInd = population[ind.sectorialindex].indiv;
			population[ind.sectorialindex].indiv = ind;//replace
			bReplaced = true;
		}
	}
	else{
		//sectorpop[ind.sectorialindex].indiv = ind;//replace
		//bReplaced = true;
		replacedInd = population[ind.sectorialindex].indiv;
		population[ind.sectorialindex].indiv = ind;
		//Pareto_HyperVolume_compare_sectorialgrid(replacedInd);// , replacedInd);
		bReplaced = true;
	}

	return bReplaced;
}

bool CHEA::update_idealpoint_from_outside()
{
	bool isUpdate = false;
	if (numprocs >1) {
		//ֻ�ж���̲�ʹ�ã���֤�����̿���ִ��
		for (int i = 0; i < partition_niche; i++)
		{
			int flag = 0;
			bool tmpUpdate = false;
			MPI_Test(&(req_idealpoint_recv[i]), &flag, &status_idealpoint);
			while (flag != 0)
			{//���յ������������̵���������

			 //cout << "[" << mpi_rank << "]("<<recv_idealpoint[i][0]<<","<<recv_idealpoint[i][1]<<","<<recv_idealpoint[i][2]<<")" << endl;
                if(!is_first)
                    tmpUpdate = update_ideal_point(recv_idealpoint[i]);
				//cout<<"["<<mpi_rank<<"]("<<recv_idealpoint[0]<<","<<recv_idealpoint[1]<<") old("<<ij[0]<<","<<ij[1]<<")now("<<idealpoint[0]<<","<<idealpoint[1]<<")"<<endl;

				//�����ȴ����ո���
				MPI_Irecv(recv_idealpoint[i], nobj, MPI_DOUBLE, neighbor_partitions[i], TAG_IDEALPOINT, MPI_COMM_WORLD, &(req_idealpoint_recv[i]));

				flag = 0;
				MPI_Test(&(req_idealpoint_recv[i]), &flag, &status_idealpoint);
			}
			isUpdate |= tmpUpdate;
		}
	}
	return isUpdate;
}


void CHEA::sendIdealpoint(double *buf, bool &isSend, int target, MPI_Request* request)
{
	if (numprocs >1)
	{
		//����̲�ʹ�ã���֤������Ҳ����ִ��
		//����������£�������һ�����̷��͸�����Ϣ
		//�������ȴ���һ������������Ϣ�������
		//if (isSend)
		//{
		//MPI_Wait(request, MPI_STATUS_IGNORE);
		//}
		//���Է��͸�����Ϣ
		//cout << "[" << mpi_rank << "] -- > " << target << "  :";
		for (int i = 0; i<nobj; i++)
		{
			//cout << idealpoint[i] << ",";
			buf[i] = IdealPoint[i];
		}
		//cout << endl;
		MPI_Isend(buf, nobj, MPI_DOUBLE, target, TAG_IDEALPOINT, MPI_COMM_WORLD, request);

		//����ָ����ϣ�������������ִ�н���������,�ߴ���߼���
		//cout<<"["<<mpi_rank<<"] send"<<endl;
		isSend = true;
	}
}

void CHEA::getAndUpdateIndivs(double* b_recv_indiv, int datasize, int source, MPI_Request *request)
{
	int flag = 0;
	//cout << "[" << mpi_rank << "] update indiv start" << endl;
	MPI_Test(request, &flag, &status_indiv);
    if (flag != 0)
	{//���յ����Խ��̵ĸ������
	 //�ӻ����л�ȡ����
	 //cout << "[" << mpi_rank << "] indiv receive" << endl;
		int m = 0;
		int updated_size = (int)(b_recv_indiv[m++] + 0.5);
		//cout << "[" << mpi_rank << "] updated size : " << updated_size << endl;
		for (int i = 0; i < updated_size; i++)
		{
			CHEAInd indiv;
			int index = (int)(b_recv_indiv[m++] + 0.5);
			for (int v = 0; v < nvar; v++)
				indiv.x_var[v] = b_recv_indiv[m++];
			for (int v = 0; v < nobj; v++)
				indiv.y_obj[v] = b_recv_indiv[m++];

            //cout << "[" << mpi_rank << "]  index : " << index << endl;//weightIndex <<"("<<lambda[weightIndex][0]<<","<<lambda[weightIndex][1]<<","<<lambda[weightIndex][2]<<")"<< endl;
            //cout << "[" << mpi_rank << "] " << indiv.x_var[0] << "," << indiv.x_var[1] << "," << indiv.x_var[2] << "," << indiv.x_var[3] << endl;
            //cout << "[" << mpi_rank << "](" << indiv.y_obj[0] << ","<<indiv.y_obj[1]<<","<<indiv.y_obj[2] << ")" << endl;

            if(!is_first)
            {
                //���������
		        update_extreme_point(indiv);

			    indiv.obj_index(IdealPoint, hyperplane_Intercept);
			    //if (indiv.sectorialindex != index)
			    //	continue;
			    //���¸���
			    CHEAInd tmp;
			    if(Pareto_HyperVolume_compare_sectorialgrid(indiv,tmp))
				    recordUpdateofIndiv[indiv.sectorialindex] = true;

            }
        }
		//cout << "[" << mpi_rank << "]  asda" << endl;
		MPI_Irecv(b_recv_indiv, datasize, MPI_DOUBLE, source, TAG_INDIV, MPI_COMM_WORLD, request);

		//cout << "[" << mpi_rank << "] after" << endl;

        //flag = 0;
        //MPI_Test(request, &flag, &status_indiv);
	}
	 //cout << "[" << mpi_rank << "] while times:" <<times << endl;
}

void CHEA::update_indiv_from_outside()
{
	if (numprocs > 1)
	{
        for (int i = 0; i <numprocs ; i++)
            if(i!=mpi_rank)
                getAndUpdateIndivs(recv_indiv[i], indivBufSize,i, &(req_indiv_recv[i]));
        //for (int i = 0; i < partition_niche; i++)
        //{
			//int neiProc = neighbor_partitions[i];
            //getAndUpdateIndivs(recv_indiv[neiProc], indivBufSize, neighbor_partitions[i], &(req_indiv_recv[neiProc]));
        //}
    }
}

void CHEA::evolution_tour_select_2()
{
	CHEAInd child2;
	bool isNeedUpdate = false;
	//bool isUpdate = false;
	int len = 0;
	set<int>::iterator it = partition_indexs[mpi_rank].begin();
	set<int>::iterator it_end = partition_indexs[mpi_rank].end();
	while (it != it_end)
	{
		int parent_index1, parent_index2;
		int b = len % (popsize / 7);
		if (b < size_subp_on_edge)
		{
			parent_index1 = subp_index_on_edge[b];
		}
		else
		{
			parent_index1 = tour_selection_hv(population);//(int)(rnd_uni(&rnd_uni_init)*popsize);//
			//parent_index1 = (int)(rnd_uni(&rnd_uni_init)*popsize);//
		}
		//parent_index2 = tour_selection_hv(sectorpop);

		//ѡ��
		while (1)
		{
			if (comm_i_mig_flag)
				parent_index2 = (int)(rnd_uni(&rnd_uni_init)*popsize); //tour_selection_hv(population);
			else
				parent_index2 = tour_selection_hv(population);//(int)(rnd_uni(&rnd_uni_init)*popsize);//
			if (parent_index2 != parent_index1)
				break;
			//cout << "[" << mpi_rank << "] index : "<<parent_index2 << endl;
		}

		//cout << "[" << mpi_rank << "] p1 = " << parent_index1 << " ; p2 = " << parent_index2 << endl;
		//����
		real_sbx_xoverA(population[parent_index1].indiv, population[parent_index2].indiv, onechild, child2);
		//����
		realmutation(onechild, 1.0 / nvar);
		//��ȡ�����Ŀ��ֵ
		onechild.obj_eval();
		//��������㡢�ο��㡢�������������Ӧ��ϵ
		//isUpdate |= update_extreme_point(onechild);
		isNeedUpdate |= update_extreme_point(onechild);

		onechild.obj_index(IdealPoint, hyperplane_Intercept);
		//Pareto_HyperVolume_compare_sectorialgrid(onechild);

		//����������
		CHEAInd tmp;
		if (Pareto_HyperVolume_compare_sectorialgrid(onechild, tmp))
			recordUpdateofIndiv[onechild.sectorialindex] = true;
		len++;
		it++;
	}
	if (isNeedUpdate&&numprocs == 1)
		update_partition();

}

//����
void CHEA::evolution_tour_select()
{
	CHEAInd child2;
	bool isUpdate = false;

	for (int i = 0; i < popsize;i++)
	{
		int parent_index1, parent_index2;


		//ѡ��
		/*
		if (rnd_uni(&rnd_uni_init) <0.381)
		{
			parent_index1 = rnd_uni(&rnd_uni_init) * popsize;
			count_selected[parent_index1]++;
			parent_index2 = rnd_uni(&rnd_uni_init) * popsize;//tour_selection_hv(sectorpop);//
			count_selected[parent_index2]++;
		}
		else
		*/{
			parent_index1 = tour_selection_hv(population);
			parent_index2 = tour_selection_hv(population);
			//parent_index1 = hv_contributions_selection(rnd_uni(&rnd_uni_init) * popsize);
			//parent_index2 = hv_contributions_selection(rnd_uni(&rnd_uni_init) * popsize);
		}

		//����
		real_sbx_xoverA(population[parent_index1].indiv, population[parent_index2].indiv, onechild, child2);
		//����
		realmutation(onechild, 1.0 / nvar);
		//��ȡ�����Ŀ��ֵ
		onechild.obj_eval();
		//��������㡢�ο��㡢�������������Ӧ��ϵ
		//
		update_extreme_point(onechild);

		onechild.obj_index(IdealPoint, hyperplane_Intercept);
		//Pareto_HyperVolume_compare_sectorialgrid(onechild);

		//����������
		CHEAInd tmp;
		//isUpdate |=
		//Pareto_HyperVolume_compare_sectorialgrid(onechild,tmp);
		if(Pareto_HyperVolume_compare_sectorialgrid(onechild,tmp))
			update_partition();
	}
	//if (isUpdate)
	//	update_partition();
}


//����
void CHEA::evolution_tour_select_2(MGraph *G)
{
	CHEAInd child2;
	bool isNeedUpdate = false;
	//bool isUpdate = false;
	int len = 0;
	set<int>::iterator it = partition_indexs[mpi_rank].begin();
	set<int>::iterator it_end = partition_indexs[mpi_rank].end();
	while (it != it_end)
	{
		int parent_index1, parent_index2;
		int b = len % (popsize/7);
		if (b < size_subp_on_edge)
		{
			parent_index1 = subp_index_on_edge[b];
		}
		else
		{
            //parent_index1 = tour_selection_hv(population);//(int)(rnd_uni(&rnd_uni_init)*popsize);//
            parent_index1 = tour_selection_Q(population,G);//(int)(rnd_uni(&rnd_uni_init)*popsize);//
            //parent_index1 = (int)(rnd_uni(&rnd_uni_init)*popsize);//
		}
		//parent_index2 = tour_selection_hv(sectorpop);

		//ѡ��
		while(1)
		{
            if(comm_i_mig_flag)
                parent_index2 =  (int)(rnd_uni(&rnd_uni_init)*popsize); //tour_selection_hv(population);
            else
                //parent_index2 = tour_selection_hv(population);//(int)(rnd_uni(&rnd_uni_init)*popsize);//
                parent_index2 = tour_selection_Q(population,G);//(int)(rnd_uni(&rnd_uni_init)*popsize);//
				//cout << "test ";
				//system("Pause");
			if(parent_index2 != parent_index1)
				break;
			//cout << "[" << mpi_rank << "] index : "<<parent_index2 << endl;
		}

		//cout << "[" << mpi_rank << "] p1 = " << parent_index1 << " ; p2 = " << parent_index2 << endl;
		//����
		real_sbx_xoverA_graph(population[parent_index1].indiv, population[parent_index2].indiv, onechild, child2);
		//����
		realmutation_graph(onechild, G);

		//��ȡ�����Ŀ��ֵ
		onechild.obj_eval_graph(G);
		//��������㡢�ο��㡢�������������Ӧ��ϵ
		//isUpdate |= update_extreme_point(onechild);
		isNeedUpdate |= update_extreme_point(onechild);

		onechild.obj_index(IdealPoint, hyperplane_Intercept);

		//Pareto_HyperVolume_compare_sectorialgrid(onechild);
        local_search_graph(onechild,G,IdealPoint,hyperplane_Intercept);
		//��ȡ�����Ŀ��ֵ
		onechild.obj_eval_graph(G);
		//��������㡢�ο��㡢�������������Ӧ��ϵ
		//isUpdate |= update_extreme_point(onechild);
		isNeedUpdate |= update_extreme_point(onechild);

		onechild.obj_index(IdealPoint, hyperplane_Intercept);

		//����������
		CHEAInd tmp;
		if (Pareto_HyperVolume_compare_sectorialgrid(onechild, tmp))
            recordUpdateofIndiv[onechild.sectorialindex] = true;
		len++;
		it++;
	}
    if (isNeedUpdate&&numprocs==1)
        update_partition();

}


void CHEA::save_front(char saveFilename[1024])
{
	std::fstream fout;
	fout.open(saveFilename, std::ios::out);
	for (int n = 0; n<population.size(); n++)
	{
		for (int k = 0; k<nobj; k++)
			fout << population[n].indiv.y_obj[k] << "  ";
		fout << "\n";
	}
	fout.close();
}

void CHEA::save_pos(char saveFilename[1024])
{
	std::fstream fout;
	fout.open(saveFilename, std::ios::out);
	for (int n = 0; n<population.size(); n++)
	{
		for (int k = 0; k<nvar; k++)
			fout << population[n].indiv.x_var[k] << "  ";
		fout << "\n";
	}
	fout.close();
}


void CHEA::operator=(const CHEA &emo)
{
	popsize = emo.popsize;
	population = emo.population;
}

void CHEA::gather_pop_y()
{
	//���͸��������Ŀ�����������ڼ���IGD��HVָ��
	int m = 0;
    set<int>::iterator it = partition_indexs[mpi_rank].begin();
    set<int>::iterator it_end = partition_indexs[mpi_rank].end();
    //cout << "[" << mpi_rank << "] partition index:";
    while (it != it_end)
    {
        int index = *it;
        //cout << index << ",";
        //for(int j =0;j<nvar;j++)
        //b_send[m++] = population[i].indiv.x_var[j];
        for (int j = 0; j<nobj; j++) {
            b_send[m++] = population[index].indiv.y_obj[j];
        }
        it++;
    }
    //for(int i=0;i<popsize;i++)
        //for (int j = 0; j<nobj; j++)
            //c_send[m++] = population[i].indiv.y_obj[j];
	//cout << endl;
    int size_send = nobj * pops_sub_self[mpi_rank];
    //int size_send = nobj * popsize;

	////cout << "["<<mpi_rank<<"]  "<<size_send<<"\t"<<m << endl;
	//�ռ�����
	int *rcounts = new int[numprocs];
	for (int i = 0; i<numprocs; i++)
        //rcounts[i] = nobj*popsize;
        rcounts[i] = nobj*pops_sub_self[i];
	int *displs = new int[numprocs];
	displs[0] = 0;
	for (int i = 1; i<numprocs; i++)
		displs[i] = displs[i - 1] + rcounts[i - 1];

	//cout << "[" << mpi_rank << "] before gatherv" << endl;
	//MPI_Gatherv(c_send, size_send, MPI_DOUBLE, c_recv, rcounts, displs, MPI_DOUBLE, root_rank, MPI_COMM_WORLD);
    MPI_Gatherv(b_send, size_send, MPI_DOUBLE, b_recv, rcounts, displs, MPI_DOUBLE, root_rank, MPI_COMM_WORLD);
	delete []rcounts;
    rcounts = NULL;
	delete []displs;
    displs = NULL;
	//cout << "[" << mpi_rank << "] gatherv end" << endl;
}

void CHEA::gather_populations()
{
	//cout << "["<<mpi_rank<<"]mark" << endl;
	//���͸��������Ŀ�����������ڼ���IGD��HVָ��
	int m = 0;
	set<int>::iterator it = partition_indexs[mpi_rank].begin();
	set<int>::iterator it_end = partition_indexs[mpi_rank].end();
	while (it != it_end)
		//for (int i = 0; i < pops_sub_self[mpi_rank]; i++)
	{
		int index = *it;
		for (int j = 0; j < nvar; j++)
			b_send[m++] = population[index].indiv.x_var[j];
		for (int j = 0; j < nobj; j++) {
			b_send[m++] = population[index].indiv.y_obj[j];
		}
		it++;
	}
	//cout <<"["<<mpi_rank<<"]"<<s_send<< "  adasd:"<<LLSize << "\t" << LLSize + pops_sub_self[mpi_rank] << endl;
	//�ռ�����
	int *rcounts = new int[numprocs];
	for (int i = 0; i<numprocs; i++)
		rcounts[i] = (nvar + nobj)*pops_sub_self[i];
	int *displs = new int[numprocs];
	displs[0] = 0;
	for (int i = 1; i<numprocs; i++)
		displs[i] = displs[i - 1] + rcounts[i - 1];
	MPI_Gatherv(b_send, s_send, MPI_DOUBLE, b_recv, rcounts, displs, MPI_DOUBLE, root_rank, MPI_COMM_WORLD);
	delete []rcounts;
	rcounts = NULL;
	delete []displs;
	displs = NULL;
	//cout << "["<<mpi_rank<<"]mark end" << endl;
}

void CHEA::calc_qualities(double &igd, double& hv)
{//ֻ��Root���̲�ʹ��
 //��������
	vector<vector<double> > pop_y;
    //pop_y.reserve(popsize);
    //pop_y.reserve(popsize*numprocs);
	int m = 0;
    for (int i = 0; i<popsize; i++)
    //for (int i = 0; i<popsize*numprocs; i++)
	{
		vector<double> y_obj;
		for (int j = 0; j<nobj; j++)
			//y_obj.push_back(c_recv[m++]);
			y_obj.push_back(b_recv[m++]);
		//cout << "[" << mpi_rank << "](" << y_obj[0] << "," << y_obj[1] << "," << y_obj[2] << ")" << endl;;
		pop_y.push_back(y_obj);
	}
	//����ȡ��ǰ��
	vector <vector<double> >  pf_y;
	population2front(pop_y, pf_y);
	//����IGD
	igd = calc_distance(pf_y);
	//����HV
	//�Ȳ���
	hv = calc_hypervolume(pf_y);

	pop_y.clear();
	pf_y.clear();
}


void CHEA::calc_qmetric(double&qmetric,int &theBestOne, MGraph *G)
{
		//vector<vector<double> > pop_x;
		//pop_x.reserve(popsize);
		//pop_y.reserve(popsize);
		//int m = 0;
		//for (int i = 0; i<popsize; i++)
		//{
			//vector<double> x_var;
			//for (int j = 0; j<nvar; j++)
				//x_var.push_back(b_recv[m++]);
			//pop_x.push_back(x_var);
			//vector<double> y_obj;
			//for (int j = 0; j<nobj; j++)
				//y_obj.push_back(b_recv[m++]);
			//pop_y.push_back(y_obj);
		//}

        qmetric=-1e+6;
        theBestOne=0;

        for(int i=0; i < popsize; i++)
        {
             //vector<vector<double> >model;
             //getModel(&b_recv[i*(nvar+nobj)],model);

             double tmp_qmetric=0.0;
             //int model_size=model.size();
             //for(int j=0; j < model_size; j++ )
             //{
                 //double q1=(getL_NRA(model[j],model[j],G)/2)/G->numEdges;
                 //int kDegree=getDegree(model[j],G);
                 //double q2=(float)(kDegree*kDegree)/(4*G->numEdges*G->numEdges);

                 //tmp_qmetric+=q1-q2;
             //}
            tmp_qmetric=-b_recv[(i+1)*(nvar+nobj)-2]-b_recv[(i+1)*(nvar+nobj)-1];
            if(tmp_qmetric>qmetric)
            {
                qmetric=tmp_qmetric;
                theBestOne=i;
            }
        }
}
double CHEA::calc_distance(vector<vector<double> >& ParetoFront)
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

double CHEA::calc_hypervolume(vector<vector<double> >& ParetoFront)
{
	wfg hv;

	double result = hv.compute(ParetoFront, hv_ref_ub);

	return result / cubeHV;
}

//���������ڱ�����ָ�������ֵ
double  CHEA::tour_selection_hv_difference(int p, vector <CHEA_SOP>  &population)
{
	int num = 0;
	int index;
	double hv_side, hv_difference = 0;

	//while (population[p].sectorialindex != population[p].indiv.sectorialindex)
	//{
		//����,��Ϊ֮ǰupdate_pop_index�Ѿ����ù�Pareto_HyperVolume_compare_sectorialgrid,���������ٱȽ�
		//p = population[p].indiv.sectorialindex;
	//}

	CHEA_SOP &subproblem = population[p];
	int size_neighborhood = subproblem.neighborhood.size();

	subproblem.indiv.obj_index(IdealPoint,hyperplane_Intercept);
	if (subproblem.indiv.sectorialindex != subproblem.sectorialindex)
		return 0;
	double hv0 = GetHyperVolume(subproblem.indiv,ReferencePoint);//

	//�����������ھ������������ֵ֮������ƽ��
	for (int i = 0; i < size_neighborhood; i++)
	{
		CHEA_SOP &nei_subp = population[subproblem.neighborhood[i]];
		if (nei_subp.sectorialindex == nei_subp.indiv.sectorialindex)
		{
			hv_side = GetHyperVolume(nei_subp.indiv, ReferencePoint);//
			hv_difference += (hv0 - hv_side) ;
			num++;
		}
	}
	if(num !=0)
		hv_difference = hv_difference / (num);

	return hv_difference;
}



//�ؼ�ָ��������
double CHEA::GetHyperVolume(const CHEAInd&  ind,vector <double> &ref_ponit_cal)
{
	double Volume = 1;

	//��׼������
	for (int j = 0; j < nobj; j++)
	{
		Volume *= (ref_ponit_cal[j] - ind.y_obj[j]);
	}

	return Volume;
}

//�ֽ������⡢������۲�����
void CHEA::gen_SubProblems(int start_obj_index, int max_value_left, vector<int> coordinate,int &index_count)
{
	if (0 == start_obj_index || 0 == max_value_left)
	{
		index_count++;

		coordinate[start_obj_index] = max_value_left;
		//}
		CHEA_SOP subproblem;
		subproblem.V_obj_o = coordinate;

		subproblem.sectorialindex = index_count;
		population.push_back(subproblem);
		int p;
		int count = 0;
		for (p = 0; p < nobj; p++)
		{
			if (coordinate[p] == 0)
			{
				subp_index_on_edge.push_back(index_count);
				break;
			}
		}
		return;
	}

	for (int i = max_value_left; i >= 0; i--)
	{
		coordinate[start_obj_index] = i;
		gen_SubProblems(start_obj_index - 1, max_value_left - i, coordinate, index_count);
	}
}



void population2front(vector <CHEA_SOP>  &mypopulation, vector <CHEAInd>  &population_front)
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
void population2front(vector <CHEA_SOP>  &mypopulation, vector <vector<double> >  &population_front)
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

void save_population(vector <CHEAInd>  &mypopulation, char saveFilename[1024])
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

//����ͼ ��ʼ��
void CHEA::init_population(MGraph *G)
{

	population[0].indiv.rnd_init(G);

	//��ʼ����ֵ
	for (int j = 0; j<nobj; j++)
	{
		TrueNadirPoint.push_back(population[0].indiv.y_obj[j]);
		IdealPoint.push_back(population[0].indiv.y_obj[j]);
		ReferencePoint.push_back(IdealPoint[j] + 1e3 * (TrueNadirPoint[j] - IdealPoint[j]));
	}
	for (int n = 1; n < popsize; n++)
	{
		population[n].indiv.rnd_init(G);
		update_extreme_point(population[n].indiv);
	}
	//�����ʼ�����Ӧ������ı��
	update_partition();
}


//ʹ��ͼ���н���
void CHEA::execute(int mg, int irun, vector<double>& qmetric, double &runtime, PTYPE ptype, MGraph *G)
{
	cout << "community detection"<<endl;
	//��ʼ������
	// first generation
	//int gen   = 1;
	//cout << "[" << mpi_rank << "] init_uniformweight" << endl;
	init_uniformweight();

	//init_neighborhood(1);
	//cout << "[" << mpi_rank << "] init_neightborhood" << endl;
	init_neighborhood();

	//if (mpi_rank == root_rank)
	//cout<< "["<<mpi_rank<<"] maxNiche : " << maxNiche << endl;

	//cout << "[" << mpi_rank << "] neighbor size : " << partition_niche << endl;

	if (numprocs > 1) {
		//cout << "[" << mpi_rank << "] something init";
		//ֻ�ж���̲�ʹ�ã���֤������Ҳ����ִ��
		//���仺��ռ�
		//���������ʹ�õĻ�����
		//cout << "[" << mpi_rank << "]  1" << endl;
		if (ptype == IDEAL_INDIV || ptype == IDEAL)
		{
			recv_idealpoint = new double*[partition_niche];
			send_idealpoint = new double*[partition_niche];
			req_idealpoint_recv = new MPI_Request[partition_niche];
			req_idealpoint_send = new MPI_Request[partition_niche];
			//isIdealpointUpdate = new bool[partition_niche];
			isBegin_Ideal = new bool[partition_niche];
		}
		//cout << "[" << mpi_rank << "]  2" << endl;
		//������½���
		//Ȩ��index + ��������+Ŀ������
		s_indiv = 1 + nvar + nobj;
		indivBufSize = 1 + s_indiv * popsize;//sizeof(int)+(sizeof(int) + sizeof(double)*(nvar+nobj)) * maxNiche;
		if (ptype == IDEAL_INDIV || ptype == INDIV)
		{
			recv_indiv = new double*[numprocs];
			send_indiv = new double*[numprocs];
			req_indiv_recv = new MPI_Request[numprocs];
			req_indiv_send = new MPI_Request[numprocs];
			isBegin_Indiv = new bool[numprocs];
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
				//������������������
				MPI_Irecv(recv_idealpoint[i], nobj, MPI_DOUBLE, neighbor_partitions[i], TAG_IDEALPOINT, MPI_COMM_WORLD, &(req_idealpoint_recv[i]));

			}
		}
		for (int i = 0; i<numprocs; i++)
		{
			if (ptype == IDEAL_INDIV || ptype == INDIV)
			{
				recv_indiv[i] = new double[indivBufSize];
				send_indiv[i] = new double[indivBufSize];
				isBegin_Indiv[i] = false;
				if (i != mpi_rank)
					MPI_Irecv(recv_indiv[i], indivBufSize, MPI_DOUBLE, i, TAG_INDIV, MPI_COMM_WORLD, &(req_indiv_recv[i]));
			}
		}
		//cout << "[" << mpi_rank << "]  4" << endl;
		MPI_Barrier(MPI_COMM_WORLD);
	}

	is_first = true;
	//cout << "[" << mpi_rank << "]  5" << endl;
	isIdealpointUpdate = false;
	recordUpdateofIndiv = new bool[population.size()];
	for (int i = 0; i < population.size(); i++)
		recordUpdateofIndiv[i] = false;

	IdealPoint.clear();
	//cout << "[" << mpi_rank << "] init population" << endl;

	init_population(G);

	//system("Pause");

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

			//rank = 0 �ռ������
			//�ռ�����
			int *rcounts = new int[numprocs];
			for (int i = 0; i < numprocs; i++) {
				rcounts[i] = nobj;
			}
			int *displs = new int[numprocs];
			displs[0] = 0;
			for (int i = 1; i < numprocs; i++)
				displs[i] = displs[i - 1] + rcounts[i - 1];

			MPI_Gatherv(send_idealpoint[0], nobj, MPI_DOUBLE, recv_utopian, rcounts, displs, MPI_DOUBLE, root_rank, MPI_COMM_WORLD);

			//ͳһ�����
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
			//�ַ���ͬ�������
			for (int i = 0; i < nobj; i++)
			{
				send_idealpoint[0][i] = IdealPoint[i];
			}
			MPI_Bcast(send_idealpoint[0], nobj, MPI_DOUBLE, root_rank, MPI_COMM_WORLD);

			//���½����������
			for (int i = 0; i < nobj; i++) {
				IdealPoint[i] = send_idealpoint[0][i];
			}
			//cout << "[" << mpi_rank << "] after sync idealpoint("<<idealpoint[0]<<","<<idealpoint[1]<<","<<idealpoint[2]<<")" << endl;
			if (mpi_rank == root_rank)
				delete[]recv_utopian;
			delete[]rcounts;
			delete[]displs;
		}
		//MPI_Barrier(MPI_COMM_WORLD);
	}

	int gen = 1;
	double qmetric_gen = 0;
    int theBestOne=0;
	char filename[1024];
	clock_t begintime;
	clock_t endtime;
	unevolvetime = 0;
	if (IS_DEBUG) {
		begintime = clock();
		//�ռ�����⼯�����ı�Ҫ����,Ŀ������
		//�ȴ���ͬ��
		if (mpi_rank == 0)
		{
			//������ʵǰ��
			//sprintf(filename, "%s/TestData/PF_Real/%s(%d).dat", WORKINGDIR, strTestInstance, nobj);
			//loadpfront(filename, ps);
		}
		if (MODE_DEBUG)
		{
			//cout << "[" << mpi_rank << "] before gather pop y" << endl;
			//gather_pop_y();
            gather_populations();
			//cout << "[" << mpi_rank << "] after" << endl;
			//Root���̴�������
			if (mpi_rank == 0)
			{
				//����⼯������IGD��HV)
                //vector<vector<double> >pop_y;
				calc_qmetric(qmetric_gen,theBestOne, G);
				qmetric.push_back(gen);  qmetric.push_back(qmetric_gen);
                 cout << "gen =  " << gen << " the best Q metric "<< qmetric_gen << "  the index of the best one " <<theBestOne<<"          NRA:"<<b_recv[(theBestOne+1)*(nvar+nobj)-2]<<" RC:"<<b_recv[(theBestOne+1)*(nvar+nobj)-1]<<endl;
			}
			//MPI_Barrier(MPI_COMM_WORLD);
		}
		//�ȴ���ͬ��
		MPI_Barrier(MPI_COMM_WORLD);
		endtime = clock();
		unevolvetime += endtime - begintime;
	}
	for (gen = 2; gen <= mg; gen++)
	{
		//cout << "the " << gen << "th generation" << endl;
		//if (gen % 100 == 0){
		//	cout << "���ǵ�" << gen << "��" << endl;
		//}
		//system("Pause");
		//����
		//ʹ��tour select �Ľ�������
		//cout << "[" << mpi_rank << "] gen = " << gen << endl;
		evolution_tour_select_2(G);

		if (ptype == IDEAL_INDIV || ptype == IDEAL)
		{
			if (numprocs > 1)
			{
				//cout <<gen<< "    [" << mpi_rank << "] update idealpoint" << endl;
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
			if (numprocs > 1 && comm_i_direct_flag&&gen % comm_i_direct == 0)
			{
				//cout << "[" << mpi_rank << "] update indiv" << endl;
				//���¸���

				for (int index = 0; index<numprocs; index++)
					if (index != mpi_rank)
					{

					int updated_size = 0;
					int m = 1;
					for (int i = 0; i < popsize; i++)
					{
						//if ((partition_indexs[mpi_rank].count(i) == 0) && (recordUpdateofIndiv[i] == true))
						if (recordUpdateofIndiv[i])
						{
							CHEAInd *_indiv = &(population[i].indiv);
							if (sop_partitions[_indiv->sectorialindex][0] == index); //note [0]
							{
								indiv_buf[m++] = i;
								for (int v = 0; v < nvar; v++)
									indiv_buf[m++] = _indiv->x_var[v];
								for (int v = 0; v < nobj; v++)
									indiv_buf[m++] = _indiv->y_obj[v];
								updated_size++;
							}
						}
					}
					if (updated_size > 0)
					{
						indiv_buf[0] = updated_size;
						int send_size = 1 + updated_size * s_indiv;


						for (int k = 0; k < 1 + updated_size*(s_indiv); k++)
							send_indiv[index][k] = indiv_buf[k];
						if (isBegin_Indiv[index])
						{
							MPI_Wait(&req_indiv_send[index], &status_indiv);
						}
						MPI_Isend(send_indiv[index], send_size, MPI_DOUBLE, index, TAG_INDIV, MPI_COMM_WORLD, &(req_indiv_send[index]));

						isBegin_Indiv[index] = true;

					}
					}
			}
			if (numprocs > 1 &&comm_i_mig_flag&& gen % comm_i_mig==0 )
			{
			int updated_size = 0;
			int m = 1;
			for (int i = 0; i < popsize;i++)
			{
			if ((partition_indexs[mpi_rank].count(i) == 0) && (recordUpdateofIndiv[i] == true))
			//if (recordUpdateofIndiv[i] == true)
			{

			indiv_buf[m++] = i;
			for (int v = 0; v < nvar; v++)
			indiv_buf[m++] = population[i].indiv.x_var[v];
			for (int v = 0; v < nobj; v++)
			indiv_buf[m++] = population[i].indiv.y_obj[v];
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

			for (int k = 0; k < 1 + updated_size*(s_indiv); k++)
			send_indiv[neiProc][k] = indiv_buf[k];
			if (isBegin_Indiv[neiProc])
			{
			MPI_Wait(&req_indiv_send[neiProc], &status_indiv);
			}
			//cout << gen<< "   [" << mpi_rank << "] -->" << neiProc << "   send indiv" << endl;
			MPI_Isend(send_indiv[neiProc], send_size, MPI_DOUBLE, neiProc, TAG_INDIV, MPI_COMM_WORLD, &(req_indiv_send[neiProc]));
			isBegin_Indiv[neiProc] = true;
			}
			}
			}
			//���ø��±�־
			if ((comm_i_mig_flag&&gen%comm_i_mig == 0) || (comm_i_direct_flag&&gen%comm_i_direct == 0))
			{
				for (int i = 0; i < popsize; i++)
					recordUpdateofIndiv[i] = false;
			}
			//cout << "[" << mpi_rank << "] update indiv end" << endl;
			update_indiv_from_outside();
		}
		//MPI_Barrier(MPI_COMM_WORLD);
		if (IS_DEBUG && gen % I_OUT == 0) {

			//��ȡ����⼯�����ı�Ҫ���ݣ�Ŀ������
			//�ȴ���ͬ��
			//update_problem_from_outside();
			//MPI_Barrier(MPI_COMM_WORLD);
			//update_problem_from_outside();
			if (MODE_DEBUG)
			{
				//�ȴ���ͬ��
				//MPI_Barrier(MPI_COMM_WORLD);
				begintime = clock();
				//cout << gen<<"   [" << mpi_rank << "] before gatherv" << endl;
				//gather_pop_y();
                gather_populations();

				//Root���̴�������
				if (mpi_rank == root_rank)
				{
					//����⼯������IGD��HV)
                vector<vector<double> >pop_y;
				calc_qmetric(qmetric_gen,theBestOne, G);
				qmetric.push_back(gen);  qmetric.push_back(qmetric_gen);
                 cout << "gen =  " << gen << " the best Q metric "<< qmetric_gen << "  the index of the best one " <<theBestOne<<"          NRA:"<<b_recv[(theBestOne+1)*(nvar+nobj)-2]<<" RC:"<<b_recv[(theBestOne+1)*(nvar+nobj)-1]<<endl;
				}
				//�ȴ���ͬ��
				MPI_Barrier(MPI_COMM_WORLD);
				endtime = clock();
				unevolvetime += endtime - begintime;
			}
			is_first = false;
		}

	}

	if (numprocs > 1)
	{
		if (ptype == IDEAL_INDIV || ptype == INDIV)
		{
			isDoneRecv = new int[numprocs];
			isDoneSend = new int[numprocs];
			b_isdone = new int[numprocs];
			req_done_recv = new MPI_Request[numprocs];
			req_done_send = new MPI_Request[numprocs];

			//�����������ɲ���
			for (int i = 0; i < numprocs; i++)
				if (i != mpi_rank)
				{
				isDoneRecv[i] = 0;
				int flag;
				MPI_Irecv(&b_isdone[i], 1, MPI_INT, i, TAG_DONE, MPI_COMM_WORLD, &req_done_recv[i]);

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
					MPI_Isend(&isDoneSend[i], 1, MPI_INT, i, TAG_DONE, MPI_COMM_WORLD, &req_done_send[i]);

				}
				}


			while (1)
			{
				for (int i = 0; i < numprocs; i++)
					if (i != mpi_rank)
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
							MPI_Isend(&isDoneSend[i], 1, MPI_INT, i, TAG_DONE, MPI_COMM_WORLD, &req_done_send[i]);

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

							getAndUpdateIndivs(recv_indiv[i], indivBufSize, i, &req_indiv_recv[i]);
						}
					}

					}
				bool isFinish = true;
				for (int i = 0; i <numprocs; i++)
					if (mpi_rank != i)
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
				//�ȴ���ͬ��
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
		}
		for (int i = 0; i<numprocs; i++)
			if (i != mpi_rank)
			{
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
    if(mpi_rank==0)
	cout <<"process ID:"<<mpi_rank<<" runtime is "<< runtime<<"s"<<endl;
	if (IS_DEBUG) {
		if (!MODE_DEBUG || (MODE_DEBUG && (gen - 1) % I_OUT != 0))
		{//���һ��������
			//��ȡ����⼯�����ı�Ҫ���ݣ�Ŀ������
			//�ȴ���ͬ��
			//MPI_Barrier(MPI_COMM_WORLD);
			//gather_pop_y();
            gather_populations();
			//Root���̴�������
			if (mpi_rank == 0)
			{
				//����⼯������IGD��HV)

				calc_qmetric(qmetric_gen,theBestOne, G);
				qmetric.push_back(gen);  qmetric.push_back(qmetric_gen);
                 cout << "gen =  " << gen << " the best Q metric "<< qmetric_gen << "  the index of the best one " <<theBestOne<<"          NRA:"<<b_recv[(theBestOne+1)*(nvar+nobj)-2]<<" RC:"<<b_recv[(theBestOne+1)*(nvar+nobj)-1]<<endl;

			}
		}
		/*
		// save the final population - X space
		//ÿ�����̶����棬���ں���չʾ
		sprintf(filename, "%s/TestData/%d/POP/POP_%s_%s(%d)_%d_%dR%d.dat",WORKINGDIR,mpi_rank,STR_MOEAD,strTestInstance,nobj, popsize, total_gen, irun);
		save_front(filename);
		sprintf(filename, "%s/TestData/%d/POF/PF_%s_%s(%d)_%d_%dR%d.dat",WORKINGDIR,mpi_rank,STR_MOEAD, strTestInstance,nobj, popsize, total_gen, irun);
		vector <CIndividual>  frontpop;
		population2front(population, frontpop);
		save_population(frontpop, filename);//�������δ浵��Ⱥ�ķ���ǰ��

		frontpop.clear();

		sprintf(filename, "%s/TestData/%d/POS/POS_%s_%s(%d)_%d_%dR%d.dat",WORKINGDIR,mpi_rank,STR_MOEAD, strTestInstance,nobj, popsize, total_gen, irun);
		save_pos(filename);
		*/
	}


	/*
	cout << "[" << mpi_rank << "] idealpoint(";
	for (int i = 0; i < nobj; i++)
	cout << IdealPoint[i] << " ";
	cout << ")" << endl;
	*/
	//������
	//�ϲ�ǰ������
	//�ȴ���ͬ��
	//	MPI_Barrier(MPI_COMM_WORLD);
	//cout << "[" << mpi_rank << "] gather" << endl;
	gather_populations();
	//cout << "[" << mpi_rank << "] gather end" << endl;
	//�ȴ���ͬ��
	//cout<<"["<<mpi_rank<<"] send data"<<endl;
	//MPI_Barrier(MPI_COMM_WORLD);

	if (mpi_rank == 0)
	{
		//cout << "receive mark" << endl;
		//cout<<"["<<mpi_rank<<"] ��������"<<endl;
		//��������
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
		vector<int> index_split;
		for (int i = 0; i<numprocs; i++)
		{
			for (int j = 0; j<pops_sub_self[i]; j++)
				index_split.push_back(i);
		}

		char filename[1024];
		pop_index.clear();
		sprintf(filename, "%s/labData/POP/POP_CAEA_%s(%d)_%d_%dR%d", WORKINGDIR, strTestInstance, nobj, popsize, total_gen, irun);
		sprintf(filename, "%s_%dnp_%s_%s",filename,numprocs,getElemTypeName(paraType),getElemTypeName(paraTopology) );
		sprintf(filename, "%s.dat",filename);
		vector <vector<double> >  pf_y;
		vector<int> pf_index;
		population2front(pop_y, pf_y, pf_index);
		//��������ǰ��
		sprintf(filename, "%s/labData/POF/PF_CAEA_%s(%d)_%d_%dR%d", WORKINGDIR,strTestInstance, nobj, popsize, total_gen, irun);
		sprintf(filename, "%s_%dnp_%s_%s",filename,numprocs,getElemTypeName(paraType),getElemTypeName(paraTopology) );
		sprintf(filename, "%s.dat",filename);
		save_population_front_split(pf_y, pf_index, index_split, filename);//�������δ浵��Ⱥ�ķ���ǰ��
		//��ӡ��Ⱥ��ÿ�������ֵ�������������ͳ�Ʒ���
		pop_y.clear();
		pf_y.clear();
		index_split.clear();
		pf_index.clear();

        sprintf(filename, "%s/labData/POS/POS_CAEA_%s(%d)_%d_%dR%d", WORKINGDIR,strTestInstance, nobj, popsize, total_gen, irun);
		sprintf(filename, "%s_%dnp_%s_%s",filename,numprocs,getElemTypeName(paraType),getElemTypeName(paraTopology) );
		sprintf(filename, "%s.dat",filename);
		std::fstream outfile;
		outfile.open(filename, std::ios::out);
        for(int i=0; i<pop_x.size();i++ )
        {
            outfile<<i;
            for(int j=0; j<pop_x[i].size();j++)
                outfile<<" "<<pop_x[i][j];
            outfile<<endl;
        }
        sprintf(filename, "%s/labData/result/result_CAEA_%s(%d)_%d_%dR%d", WORKINGDIR,strTestInstance, nobj, popsize, total_gen, irun);
		sprintf(filename, "%s_%dnp_%s_%s",filename,numprocs,getElemTypeName(paraType),getElemTypeName(paraTopology) );
		sprintf(filename, "%s.dat",filename);
		std::fstream fout;
		fout.open(filename, std::ios::out);
		for (int input_time = 0; input_time < popsize; input_time++){
			fout << "community-detection-result: "<<input_time<<"\n";
			vector<double > result = pop_x[0+input_time];
			vector<vector<double> > model;
			getModel(result,model);
			for (int i = 0; i < model.size(); i++){
                int j=0;
                if(model[i][0]==model[i][1]&&model[i][0]==0)j=1;
				for ( ;j < model[i].size(); j++){
					if (j == model[i].size()-1)
					{
						fout << model[i][j] ;
					}
					else{
						fout << model[i][j] << " ";
					}
				}
				fout << "\n";
			}
		}

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
		}
		for (int i = 0; i<numprocs; i++)
		{
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

		sop_partitions.clear();
		partition_indexs.clear();
		partition_overlap_indexs.clear();
		neighbor_partitions_set.clear();
		neighbor_partitions.clear();


	}
}
