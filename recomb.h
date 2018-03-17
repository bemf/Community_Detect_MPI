#ifndef _RECOMBINATION_H_
#define _RECOMBINATION_H_

#ifndef __GLOBAL_H_
#include "global.h"
#endif
#ifndef TESTINSTANCE_H
#include "TestFunction/testInstance.h"
#endif
//#include "individual.h"

/* Routine for real polynomial mutation of an T */

template <class T>
void realmutation(T &ind, double rate)
{
	double rnd, delta1, delta2, mut_pow, deltaq;
	double y, yl, yu, val, xy;
	double eta_m = etam;

	int    id_rnd = int(rnd_uni(&rnd_uni_init)*nvar);

	for (int j=0; j<nvar; j++)
	{
		if (rnd_uni(&rnd_uni_init)<=rate)
		{
			y  = ind.x_var[j];
			yl = lowBound[j];
			yu = uppBound[j];
			delta1 = (y-yl)/(yu-yl);
			delta2 = (yu-y)/(yu-yl);
			rnd = rnd_uni(&rnd_uni_init);
			mut_pow = 1.0/(eta_m+1.0);
			if (rnd <= 0.5)
			{
				xy = 1.0-delta1;
				val = 2.0*rnd+(1.0-2.0*rnd)*(pow(xy,(eta_m+1.0)));
				deltaq =  pow(val,mut_pow) - 1.0;
			}
			else
			{
				xy = 1.0-delta2;
				val = 2.0*(1.0-rnd)+2.0*(rnd-0.5)*(pow(xy,(eta_m+1.0)));
				deltaq = 1.0 - (pow(val,mut_pow));
			}
			y = y + deltaq*(yu-yl);
			if (y<yl)
				//y = yl;
				y = yl + ( (yl-y)/(yu-yl) - floor( (yl-y)/(yu-yl) )  ) * (yu-yl);
			if (y>yu)
				//y = yu;
				y = yu - ( (y-yu)/(yu-yl) - floor( (y-yu)/(yu-yl) )  ) * (yu-yl);
			ind.x_var[j] = y;
		}
	}
	return;
}

template <class T>
void diff_evo_xoverB(T &ind0, T &ind1, T &ind2, T &child, double rate, double CR)
{
	int idx_rnd = int(rnd_uni(&rnd_uni_init)*nvar);

	for(int n=0;n<nvar;n++)
	{
	  /*Selected Two Parents*/

	  // strategy one
	  // child.x_var[n] = ind0.x_var[n] + rate*(ind2.x_var[n] - ind1.x_var[n]);

	  ///*
	  // strategy two
	  double rnd1 = rnd_uni(&rnd_uni_init);

	  if(rnd1<CR||n==idx_rnd)
		  child.x_var[n] = ind0.x_var[n] + rate*(ind2.x_var[n] - ind1.x_var[n]);
	  else
		  child.x_var[n] = ind0.x_var[n];
	  //*/


	  // handle the boundary voilation
	  if(child.x_var[n]<lowBound[n]){
		  //double rnd = rnd_uni(&rnd_uni_init);
		  //child.x_var[n] = lowBound[n] + rnd*(ind0.x_var[n] - lowBound[n]);
		  //child.x_var[n] = 2 * lowBound[n] - child.x_var[n];
		  child.x_var[n] = lowBound[n] + ( (lowBound[n]-child.x_var[n])/(uppBound[n]-lowBound[n]) - floor( (lowBound[n]-child.x_var[n])/(uppBound[n]-lowBound[n]) )  ) * (uppBound[n]-lowBound[n]);
	  }
	  if(child.x_var[n]>uppBound[n]){
		  //double rnd = rnd_uni(&rnd_uni_init);
		  //child.x_var[n] = uppBound[n] - rnd*(uppBound[n] - ind0.x_var[n]);
		  //child.x_var[n] = 2 * uppBound[n] - child.x_var[n];
		  child.x_var[n] = uppBound[n] - ( (child.x_var[n]-uppBound[n])/(uppBound[n]-lowBound[n]) - floor( (child.x_var[n]-uppBound[n])/(uppBound[n]-lowBound[n]) )  ) * (uppBound[n]-lowBound[n]);
	  }

	  //if(child.x_var[n]<lowBound) child.x_var[n] = lowBound;
	  //if(child.x_var[n]>uppBound) child.x_var[n] = uppBound;
	}
}


template <class T>
void realmutation_graph(T &ind, MGraph *G){

    int rand_i = 0;
    //rand_i = int(rnd_uni(&rnd_uni_init) * nvar);
    if(rnd_uni(&rnd_uni_init)<=graph_mutation_p)
    {
        //for(int j=0;j<10;j++)
        {
        rand_i = int(rnd_uni(&rnd_uni_init) * nvar);
        int Num = 0;
        int rand_temp = int(rnd_uni(&rnd_uni_init) *G->degree[rand_i]) ;//从1到Num生成随机数
        ind.x_var[rand_i] = G->Adjacency[rand_i][rand_temp];




        }
    }

	return;
}

template <class T>
void local_search_graph(T &ind, MGraph *G,const vector<double>&idealPoint,const int &hyperplane_Intercept, int MI=1)
{
    vector<int>V_obj(nobj,0);
    ind.cal_V_obj(idealPoint,hyperplane_Intercept,V_obj);
    int count =0;
    bool ischanged=true;
    while(count<MI&&ischanged==true)
   // while(ischanged==true)
    {
        ischanged=false;
        count++;

        vector<int> index;
        index.reserve(nvar);
        for(int i=0;i<nvar;i++)index.push_back(i);

        int n=nvar;
        vector<vector<int> >model;
        vector<int>community_label;
        getModel(ind.x_var,model,community_label);
        for(int i=0;i<nvar;i++)
        {
            int rnd_i=int(rnd_uni(&rnd_uni_init)*n);
            rnd_i=index[rnd_i];

            double max_delta_Q=0.0;
            int theBestOne=-1;
            int i_label=community_label[rnd_i];
            int j_label;
            int ASize=G->degree[rnd_i];
            int l_vC=0,Kc=0,Kv=ASize;

            int i_size=model[i_label].size();
            for(int k=0;k<i_size;k++)
                if(model[i_label][k]!=rnd_i)
                    Kc+=G->degree[model[i_label][k]];
            for(int k=0;k<i_size;k++)
                if(model[i_label][k]!=rnd_i)
                    l_vC+=G->arc[model[i_label][k]][rnd_i];
            //max_delta_Q=(double)(l_vC)/G->numEdges+(double)(-Kc*Kv)/(2*G->numEdges*G->numEdges);
            max_delta_Q=(double)(l_vC)/G->numEdges*V_obj[0]+(double)(-Kc*Kv)/(2*G->numEdges*G->numEdges)*V_obj[1];

            for(int j=0;j<ASize;j++)
            {
                int node=G->Adjacency[rnd_i][j];
                j_label=community_label[node];
                if(i_label!=j_label)
                {
                    int l_vB,Kb;
                    l_vB=Kb=0;
                    int j_size=model[j_label].size();
                    for(int k=0;k<j_size;k++)Kb+=G->degree[model[j_label][k]];
                    for(int k=0;k<j_size;k++)l_vB+=G->arc[model[j_label][k]][rnd_i];




                    double tmp_delta_Q=(double)(l_vB)/G->numEdges*V_obj[0]+(double)(-Kb*Kv)/(2*G->numEdges*G->numEdges)*V_obj[1];
                    //double tmp_delta_Q=(double)(l_vB)/G->numEdges+(double)(-Kb*Kv)/(2*G->numEdges*G->numEdges);
                    if(max_delta_Q<tmp_delta_Q)
                    {
                         max_delta_Q=tmp_delta_Q;
                         theBestOne=node;
                    }
                }
            }
            if(theBestOne!=-1&&community_label[theBestOne]!=i_label)
            {
                 ind.x_var[rnd_i]=theBestOne;

                 //int size=model[i_label].size();
                 for(int j=0;j<i_size;j++)
                     if(model[i_label][j]==rnd_i)model[i_label].erase(model[i_label].begin()+j);
                 model[community_label[theBestOne]].push_back(rnd_i);

                 community_label[rnd_i]=community_label[theBestOne];
                 ischanged=true;

            }

            int tmp=index[rnd_i];
            index[rnd_i]=index[n-1];
            index[n-1]=tmp;
            n--;
        }

    }
    //double local_search_rate=1;
    //double max_qmetric=-1e+6;
    //int theBestOne=0;
    //if(rnd_uni(&rnd_uni_init)<=local_search_rate)
    //{
        ////cout<<"local_search"<<endl;
        //int size=G->degree[rand_i];
        //for (int i = 0; i < size; i++)
        //{

                //vector<double>x=ind.x_var;
                //x[rand_i]=G->Adjacency[rand_i][i];

                //vector<vector<double> >model;
                //getModel(x,model);

                //double tmp_qmetric=0.0;
                //int model_size=model.size();
                //for(int j=0; j < model_size; j++ )
                //{
                    ////if(model[j][0]==model[j][1]&&model[j][1]==0)model[j].erase(model[j].begin());
                    //double q1=(getL_NRA(model[j],model[j],G)/2)/G->numEdges;
                    ////double q2=(float)(getDegree(model[j],G)*getDegree(model[j],G))/(4*G->numEdges*G->numEdges);
                    //int kDegree=getDegree(model[j],G);
                    //double q2=(float)(kDegree*kDegree)/(4*G->numEdges*G->numEdges);
                    //tmp_qmetric+=q1-q2;
                 //}
                //if(tmp_qmetric>max_qmetric)
                //{
                    //max_qmetric=tmp_qmetric;
                    //theBestOne=G->Adjacency[rand_i][i];
               //}

        //}
        //ind.x_var[rand_i]=theBestOne;
    //}
    return ;
}
/* Routine for real variable SBX crossover */
template <class T>
void real_sbx_xoverA_graph(T &parent1, T &parent2, T &child1, T &child2){
	int rand_i = 1;
	int rand_j = 0;
	rand_i = int(rnd_uni(&rnd_uni_init) * nvar);
	rand_j = int(rnd_uni(&rnd_uni_init) * (nvar - rand_i)) + rand_i;
	double y1, y2;

    if(rnd_uni(&rnd_uni_init)<=graph_xcross_p)
    {
        for(int i=0; i< rand_i; i++)
        {
    		child1.x_var[i] = parent1.x_var[i];
    		child2.x_var[i] = parent2.x_var[i];
        }
    	for (int i = rand_i; i <= rand_j; i++){
    		y1 = parent1.x_var[i];
    		y2 = parent2.x_var[i];
    		child1.x_var[i] = y2;
    		child2.x_var[i] = y1;
    	}
        for(int i=rand_j+1; i< nvar; i++)
        {
    		child1.x_var[i] = parent1.x_var[i];
    		child2.x_var[i] = parent2.x_var[i];
        }
    }
    else
    {
    	for (int i=0; i<nvar; i++)
    	{
    		child1.x_var[i] = parent1.x_var[i];
    		child2.x_var[i] = parent2.x_var[i];
    	}
    }
	return;
}
template <class T>
void real_sbx_xoverA(T &parent1, T &parent2, T &child1, T &child2)
{
	double rand;
	double y1, y2, yl, yu;
	double c1, c2;
	double alpha, beta, betaq;
	double eta_c = etax;
	if (rnd_uni(&rnd_uni_init) <= 1.0)
	{
		for (int i=0; i<nvar; i++)
		{
			if (rnd_uni(&rnd_uni_init)<=0.5 )
			{
				if (fabs(parent1.x_var[i]-parent2.x_var[i]) > EPS)
				{
					if (parent1.x_var[i] < parent2.x_var[i])
					{
						y1 = parent1.x_var[i];
						y2 = parent2.x_var[i];
					}
					else
					{
						y1 = parent2.x_var[i];
						y2 = parent1.x_var[i];
					}
					yl = lowBound[i];
					yu = uppBound[i];
					rand = rnd_uni(&rnd_uni_init);
					beta = 1.0 + (2.0*(y1-yl)/(y2-y1));
					alpha = 2.0 - pow(beta,-(eta_c+1.0));
					if (rand <= (1.0/alpha))
					{
						betaq = pow ((rand*alpha),(1.0/(eta_c+1.0)));
					}
					else
					{
						betaq = pow ((1.0/(2.0 - rand*alpha)),(1.0/(eta_c+1.0)));
					}
					c1 = 0.5*((y1+y2)-betaq*(y2-y1));
					beta = 1.0 + (2.0*(yu-y2)/(y2-y1));
					alpha = 2.0 - pow(beta,-(eta_c+1.0));
					if (rand <= (1.0/alpha))
					{
						betaq = pow ((rand*alpha),(1.0/(eta_c+1.0)));
					}
					else
					{
						betaq = pow ((1.0/(2.0 - rand*alpha)),(1.0/(eta_c+1.0)));
					}
					c2 = 0.5*((y1+y2)+betaq*(y2-y1));
					if (c1<yl)
						//c1=yl;
						c1 = yl + ( (yl-c1)/(yu-yl) - floor( (yl-c1)/(yu-yl) )  ) * (yu-yl);
					if (c2<yl)
						//c2=yl;
						c2 = yl + ( (yl-c2)/(yu-yl) - floor( (yl-c2)/(yu-yl) )  ) * (yu-yl);
					if (c1>yu)
						//c1=yu;
						c1 = yu - ( (c1-yu)/(yu-yl) - floor( (c1-yu)/(yu-yl) )  ) * (yu-yl);
					if (c2>yu)
						//c2=yu;
						c2 = yu - ( (c2-yu)/(yu-yl) - floor( (c2-yu)/(yu-yl) )  ) * (yu-yl);
					if (rnd_uni(&rnd_uni_init)<=0.5)
					{
						child1.x_var[i] = c2;
						child2.x_var[i] = c1;
					}
					else
					{
						child1.x_var[i] = c1;
						child2.x_var[i] = c2;
					}
				}
				else
				{
					child1.x_var[i] = parent1.x_var[i];
					child2.x_var[i] = parent2.x_var[i];
				}
			}
			else
			{
				child1.x_var[i] = parent1.x_var[i];
				child2.x_var[i] = parent2.x_var[i];
			}
		}
	}
	else
	{
		for (int i=0; i<nvar; i++)
		{
			child1.x_var[i] = parent1.x_var[i];
			child2.x_var[i] = parent2.x_var[i];
		}
	}
	return;
}


/* Routine for real variable SBX crossover */
//template <class T>
//void realbinarycrossover (T &parent1, T &parent2, T& child1, T& child2)
//{
//	double rand;
//	double y1, y2, yl, yu;
//	double c1, c2;
//	double alpha, beta, betaq;
//	double eta_c = etax;
//	if (rnd_uni(&rnd_uni_init) <= 1.0)
//	{
//		for (int i=0; i<nvar; i++)
//		{
//			if (rnd_uni(&rnd_uni_init)<=0.5 )
//			{
//				if (fabs(parent1.x_var[i]-parent2.x_var[i]) > EPS)
//				{
//					if (parent1.x_var[i] < parent2.x_var[i])
//					{
//						y1 = parent1.x_var[i];
//						y2 = parent2.x_var[i];
//					}
//					else
//					{
//						y1 = parent2.x_var[i];
//						y2 = parent1.x_var[i];
//					}
//					yl = lowBound[i];
//					yu = uppBound[i];
//					rand = rnd_uni(&rnd_uni_init);
//					beta = 1.0 + (2.0*(y1-yl)/(y2-y1));
//					alpha = 2.0 - pow(beta,-(eta_c+1.0));
//					if (rand <= (1.0/alpha))
//					{
//						betaq = pow ((rand*alpha),(1.0/(eta_c+1.0)));
//					}
//					else
//					{
//						betaq = pow ((1.0/(2.0 - rand*alpha)),(1.0/(eta_c+1.0)));
//					}
//					c1 = 0.5*((y1+y2)-betaq*(y2-y1));
//					beta = 1.0 + (2.0*(yu-y2)/(y2-y1));
//					alpha = 2.0 - pow(beta,-(eta_c+1.0));
//					if (rand <= (1.0/alpha))
//					{
//						betaq = pow ((rand*alpha),(1.0/(eta_c+1.0)));
//					}
//					else
//					{
//						betaq = pow ((1.0/(2.0 - rand*alpha)),(1.0/(eta_c+1.0)));
//					}
//					c2 = 0.5*((y1+y2)+betaq*(y2-y1));
//					//if (c1<yl)
//					//	c1=yl;
//					//if (c2<yl)
//					//	c2=yl;
//					//if (c1>yu)
//					//	c1=yu;
//					//if (c2>yu)
//					//	c2=yu;
//					if (c1<yl)
//						c1 = yl + ( (yl-c1)/(yu-yl) - floor( (yl-c1)/(yu-yl) )  ) * (yu-yl);
//					if (c2<yl)
//						c2 = yl + ( (yl-c2)/(yu-yl) - floor( (yl-c2)/(yu-yl) )  ) * (yu-yl);
//					if (c1>yu)
//						c1 = yu - ( (c1-yu)/(yu-yl) - floor( (c1-yu)/(yu-yl) )  ) * (yu-yl);
//					if (c2>yu)
//						c2 = yu - ( (c2-yu)/(yu-yl) - floor( (c2-yu)/(yu-yl) )  ) * (yu-yl);
//					if (rnd_uni(&rnd_uni_init)<=0.5)
//					{
//						child1.x_var[i] = c2;
//						child2.x_var[i] = c1;
//					}
//					else
//					{
//						child1.x_var[i] = c1;
//						child2.x_var[i] = c2;
//					}
//				}
//				else
//				{
//					child1.x_var[i] = parent1.x_var[i];
//					child2.x_var[i] = parent2.x_var[i];
//				}
//			}
//			else
//			{
//				child1.x_var[i] = parent1.x_var[i];
//				child2.x_var[i] = parent2.x_var[i];
//			}
//		}
//	}
//	else
//	{
//		for (int i=0; i<nvar; i++)
//		{
//			child1.x_var[i] = parent1.x_var[i];
//			child2.x_var[i] = parent2.x_var[i];
//		}
//	}
//	return;
//}

template <class T>
void diff_evo_xover(T ind0, T ind1, T ind2, T ind3, T& new_ind, double rate)
{

	// Check Whether the cross-over to be performed
	/*Loop over no of variables*/
	int idx_rnd = int(rnd_uni(&rnd_uni_init)*nvar);

	for (int n = 0; n<nvar; n++){
		/*Selected Two Parents*/

		double rnd = rnd_uni(&rnd_uni_init);
		if (rnd<1 || n == idx_rnd)
			new_ind.x_var[n] = ind1.x_var[n] + rate*(ind2.x_var[n] - ind3.x_var[n]);
		else
			new_ind.x_var[n] = ind0.x_var[n];

		if (new_ind.x_var[n]<lowBound[n]) new_ind.x_var[n] = lowBound[n];
		if (new_ind.x_var[n]>uppBound[n]) new_ind.x_var[n] = uppBound[n];
	}
}




template <class T>
void diff_evo_xover2(T ind0, T ind1, T ind2, T& new_ind,double rate)
{

	int idx_rnd = int(rnd_uni(&rnd_uni_init)*nvar);

	for (int n = 0; n<nvar; n++)
	{
		/*Selected Two Parents*/

		new_ind.x_var[n] = ind0.x_var[n] + rate*(ind2.x_var[n] - ind1.x_var[n]);

		/*
		double rnd1 = rnd_uni(&rnd_uni_init);
		if(rnd1<0.9||n==idx_rnd)
		//if(rnd1<1.0||n==idx_rnd)
		new_ind.x_var[n] = ind0.x_var[n] + rate*(ind2.x_var[n] - ind1.x_var[n]);
		else
		new_ind.x_var[n] = ind0.x_var[n];
		//*/

		if (new_ind.x_var[n]<lowBound[n]){
			double rnd = rnd_uni(&rnd_uni_init);
			new_ind.x_var[n] = lowBound[n] + rnd*(ind0.x_var[n] - lowBound[n]);
		}
		if (new_ind.x_var[n]>uppBound[n]){
			double rnd = rnd_uni(&rnd_uni_init);
			new_ind.x_var[n] = uppBound[n] - rnd*(uppBound[n] - ind0.x_var[n]);
		}

		//if(new_ind.x_var[n]<lowBound) new_ind.x_var[n] = lowBound;
		//if(new_ind.x_var[n]>uppBound) new_ind.x_var[n] = uppBound;
	}
}

template <class T>
void diff_evo_xover3(T ind0, T ind1, T ind2, vector<double> &xdiff, T& new_ind, double rate)
{
	double rnd = rnd_uni(&rnd_uni_init), rnd2 = rnd_uni(&rnd_uni_init);
	for (int n = 0; n<nvar; n++)
	{
		/*Selected Two Parents*/

		if (rnd<1)
			new_ind.x_var[n] = ind0.x_var[n] + rate*(ind2.x_var[n] - ind1.x_var[n]);
		else
			new_ind.x_var[n] = ind0.x_var[n] + rnd2*xdiff[n];

		if (new_ind.x_var[n]<lowBound[n]) new_ind.x_var[n] = lowBound[n];
		if (new_ind.x_var[n]>uppBound[n]) new_ind.x_var[n] = uppBound[n];
	}
}

#endif
