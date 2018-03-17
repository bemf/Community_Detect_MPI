#include "testInstance.h"


void initTestInstance(){
    if(!strcmp(strTestInstance,"ZDT1")|| !strcmp(strTestInstance,"ZDT2") || !strcmp(strTestInstance,"ZDT3") || !strcmp(strTestInstance,"ZDT4") || !strcmp(strTestInstance,"ZDT6") || !strcmp(strTestInstance,"OKA-1") || !strcmp(strTestInstance,"OKA-2"))
	{
		lowBound.clear();
		int i=0;
		for(;i<nvar;i++)
            lowBound.push_back(0);

		uppBound.clear();
		for(i=0;i<nvar;i++)
            uppBound.push_back(1);

		hv_ref_ub.resize(2,11);
		if(!strcmp(strTestInstance,"ZDT3"))
			hv_ref_lb.resize(2,-0.8);
		else
			hv_ref_lb.resize(2,0);
	}

	else if(!strcmp(strTestInstance,"DTLZ1") || !strcmp(strTestInstance,"DTLZ2") ||!strcmp(strTestInstance,"DTLZ3") || !strcmp(strTestInstance,"DTLZ4") || !strcmp(strTestInstance,"DTLZ5") || !strcmp(strTestInstance,"DTLZ6") || !strcmp(strTestInstance,"DTLZ7") || !strcmp(strTestInstance,"Convex_DTLZ2") || !strcmp(strTestInstance,"I_DTLZ1") )
	{
		lowBound.resize(nvar);
		uppBound.resize(nvar);
		DTLZ::DTLZ_Vars_Allocate(&(*(lowBound.begin())),&(*(uppBound.begin())),nvar);

		hv_ref_ub.resize(nobj,11);
		hv_ref_lb.resize(nobj,0);
	}

	else if (!strcmp(strTestInstance, "UF1") || !strcmp(strTestInstance, "UF2") || !strcmp(strTestInstance, "UF3") || !strcmp(strTestInstance, "UF4") || !strcmp(strTestInstance, "UF5") || !strcmp(strTestInstance, "UF6") || !strcmp(strTestInstance, "UF7") || !strcmp(strTestInstance, "UF8") || !strcmp(strTestInstance, "UF9") || !strcmp(strTestInstance, "UF10") )
    {
		lowBound.resize(nvar);
		uppBound.resize(nvar);

		CEC09::CEC09_Vars_Allocate(&(*(lowBound.begin())),&(*(uppBound.begin())),nvar,CEC09::getCEC09Type(strTestInstance));
		hv_ref_ub.resize(nobj,11);
		hv_ref_lb.resize(2,0);
    }

	cubeHV = 1.0;
	for(int i = 0 ;  i < nobj;i++)
	{
		cubeHV *= (hv_ref_ub[i] - hv_ref_lb[i]);
	}

}
void  getModel(double* x,vector<vector<double> >&model)
{
    model.reserve(nvar);
    vector<vector<int> >Adjacency;
    Adjacency.resize(nvar);
    vector<int>visited(nvar,0);
    vector<int>stack;
    stack.reserve(nvar);

    for(int i=0;i<nvar;i++)
        Adjacency[i].reserve(nvar);

    for(int i=0;i<nvar;i++)
        if(x[i]!=i)
        {
         Adjacency[i].push_back(x[i]);
         Adjacency[x[i]].push_back(i);
        }

    for(int i=0;i<nvar;i++)
    {
        if(visited[i]==0)
        {
            stack.push_back(i);
            vector<double>community;
            community.reserve(nvar);
            community.push_back(i);
            visited[i]=1;
            while(!stack.empty())
            {
                int top=stack.back();
                stack.pop_back();
                ////has been visited

                int adjSize=Adjacency[top].size();
                for(int j=0;j<adjSize;j++)
                {
                    int target=Adjacency[top][j];
                     if(visited[target]==0)
                     {
                        stack.push_back(target);
                        visited[target]=1;
                        community.push_back(target);
                     }
                }

            }
            model.push_back(community);
        }
    }
}
void  getModel(vector<double> &x,vector<vector<int> >&model,vector<int>&label)
{
    model.reserve(nvar);
    vector<vector<int> >Adjacency;
    Adjacency.resize(nvar);
    vector<int>visited(nvar,0);
    vector<int>stack;
    stack.reserve(nvar);

    label.resize(nvar);
    int label_i=0;

    for(int i=0;i<nvar;i++)
        Adjacency[i].reserve(nvar);

    for(int i=0;i<nvar;i++)
        if(x[i]!=i)
        {
         Adjacency[i].push_back(x[i]);
         Adjacency[x[i]].push_back(i);
        }

    for(int i=0;i<nvar;i++)
    {
        if(visited[i]==0)
        {
            stack.push_back(i);
            vector<int>community;
            community.reserve(nvar);
            community.push_back(i);
            label[i]=label_i;
            visited[i]=1;
            while(!stack.empty())
            {
                int top=stack.back();
                stack.pop_back();
                ////has been visited

                int adjSize=Adjacency[top].size();
                for(int j=0;j<adjSize;j++)
                {
                    int target=Adjacency[top][j];
                     if(visited[target]==0)
                     {
                        stack.push_back(target);
                        visited[target]=1;
                        community.push_back(target);
                        label[target]=label_i;
                     }
                }

            }
            model.push_back(community);
            label_i++;
        }
    }


}
void  getModel(vector<double> &x,vector<vector<double> >&model)
{
    model.reserve(nvar);
    vector<vector<int> >Adjacency;
    Adjacency.resize(nvar);
    vector<int>visited(nvar,0);
    vector<int>stack;
    stack.reserve(nvar);

    for(int i=0;i<nvar;i++)
        Adjacency[i].reserve(nvar);

    for(int i=0;i<nvar;i++)
        if(x[i]!=i)
        {
         Adjacency[i].push_back(x[i]);
         Adjacency[x[i]].push_back(i);
        }

    for(int i=0;i<nvar;i++)
    {
        if(visited[i]==0)
        {
            stack.push_back(i);
            vector<double>community;
            community.reserve(nvar);
            community.push_back(i);
            visited[i]=1;
            while(!stack.empty())
            {
                int top=stack.back();
                stack.pop_back();
                ////has been visited

                int adjSize=Adjacency[top].size();
                for(int j=0;j<adjSize;j++)
                {
                    int target=Adjacency[top][j];
                     if(visited[target]==0)
                     {
                        stack.push_back(target);
                        visited[target]=1;
                        community.push_back(target);
                     }
                }

            }
            model.push_back(community);
        }
    }


}
//将现有的某个个体的值转化成种群
vector<vector<double> > getModel(vector<double> x){
	vector<vector<double> > model;
	vector<double>::iterator iter_i;
	vector<double>::iterator iter_x;
	vector<double> temp_x;
	for (int i = 0; i < x.size(); i++){
		if (i == 0){
			vector<double> t;
			/*t.push_back(i);
			t.push_back(x[i]);*/
			addtemp(t, i);
			addtemp(t, x[i]);
			model.push_back(temp_x);
		}
		iter_i = find(temp_x.begin(), temp_x.end(), i);
		iter_x = find(temp_x.begin(), temp_x.end(), x[i]);
		if (iter_i == temp_x.end() && iter_x == temp_x.end()){
			vector<double> t;
			/*t.push_back(i);
			t.push_back(x[i]);*/
			addtemp(t, i);
			addtemp(t, x[i]);
			model.push_back(t);
		}
		else{
			if (iter_i != temp_x.end() && iter_x == temp_x.end())
			{
				int index = findmodel(model, i);
				insertmodel(model, x[i], index);
			}
			if (iter_x != temp_x.end() && iter_i == temp_x.end())
			{
				int index = findmodel(model, x[i]);
				insertmodel(model, i, index);
			}
			int model_index1 = findmodel(model, i);
			int model_index2 = findmodel(model, x[i]);
			if (model_index1 != model_index2) {
				for (int n = 0; n < model[model_index2].size();n++) {
					insertmodel(model,model[model_index2][n],model_index1);
				}
				for (int n = 0; n < model[model_index2].size(); n++) {
					model[model_index2][n] = -1;
				}
			}
		}
		addtemp(temp_x, i);
		addtemp(temp_x, x[i]);
	}
	//cout << "______________________________________________________" << endl;
	//for (int i = 0; i < model.size(); i++) {
	//	for (int j = 0; j < model[i].size(); j++) {
	//		cout << model[i][j] << " ";
	//	}
	//	cout << endl;
	//}
	//cout << "______________________________________________________" << endl;
	vector<vector<double> > dealmodel;
	vector<double> temp_add;
	for (int i = 0; i < model.size(); i++) {
		for (int j = 0; j < model[i].size(); j++) {
			if (model[i][j]!=-1)
			{
				temp_add.push_back(model[i][j]);
			}
		}
		if (temp_add.size() != 0) {
			dealmodel.push_back(temp_add);
		}
		temp_add.clear();
	}
	int num = 0;
	for (int i = 0; i < dealmodel.size(); i++) {
		for (int j = 0; j < dealmodel[i].size(); j++) {
			num++;
		}
	}
	//if (num != 34) {
		//cout << "出错了" << num<<endl;
		//cout << "x:____________________________________" << endl;
		//for (int i = 0; i < x.size(); i++) {
			//cout << x[i] << " ";
		//}
		//cout << endl;
		//cout << "x:____________________________________" << endl;
		//for (int i = 0; i < dealmodel.size(); i++) {
			//for (int j = 0; j < dealmodel[i].size(); j++) {
				//cout << dealmodel[i][j] << " ";
			//}
			//cout << endl;
		//}
		//system("Pause");
	//}
	return dealmodel;
}
void insertmodel(vector<vector<double> >&model, int p, int index) {
	bool can = true;
	for (int j = 0; j < model[index].size(); j++)
	{
		if (model[index][j] == p) {
			can = false;
		}
	}
	if (can) {
		model[index].push_back(p);
	}
}
void addtemp(vector<double>&temp,int p) {
	bool can = true;
	for (int j = 0; j < temp.size(); j++) {
		if (temp[j] == p) {
			can = false;
		}
	}
	if (can) {
		temp.push_back(p);
	}
}
int findmodel(vector<vector<double> >model, int x) {
	int index = -1;
	for (int i = 0; i < model.size(); i++) {
		for (int j = 0; j < model[i].size(); j++) {
			if (x == model[i][j]) {
				return i;
			}
		}
	}
	return index;
}
//社区检测第一个目标函数
double getL_NRA(vector<double> v1, vector<double> v2, MGraph *G){
	double l = 0;
    int v1_size=v1.size();
	for (int i = 0; i <v1_size; i++)
	{
		int g_i = v1[i];
        int v2_size=v2.size();
		for (int j = 0; j < v2_size; j++)
		{
			int g_j = v2[j];
			l = l + G->arc[g_i][g_j];
		}
	}
	return l;
}
//社区检测第二个目标函数
double getL_RC(vector<double> v1, MGraph *G){
	double l = 0;
	vector<double> v2;
	vector<double>::iterator iter;
	for (int j = 0; j < G->numVertexes; j++)
	{
		iter = find(v1.begin(), v1.end(), j);
		if (iter == v1.end()){
			v2.push_back(j);
		}
	}
	//for (int i = 0; i < v1.size(); i++){
	//	cout << v1[i]<<",";
	//}
	//cout << endl;
	//for (int i = 0; i < v2.size(); i++){
	//	cout << v2[i] << ",";
	//}

	for (int i = 0; i < v1.size(); i++)
	{
		int g_i = v1[i];
		for (int j = 0; j <v2.size(); j++)
		{
			int g_j = v2[j];
			l = l + G->arc[g_i][g_j];
		}
	}
	return l;
}

int getDegree(vector<double> &m,MGraph *G)
{
    int d=0;
    int m_size=m.size();
    for(int i=0; i<m_size; i++)
    {
        d+=G->degree[(int)m[i]];
    }
    return d;
}

//社区检测中计算两个目标函数的值
void objectives_graph1(vector<double> &x_var, vector <double> &y_obj,MGraph *G){
	vector<vector<double> >model ;
	getModel(x_var,model);
	//y_graph[0],NRA
	double NRA = 0;
	for (int i = 0; i < model.size(); i++){
        if(model[i][0]==model[i][1]&&model[i][1]==0)model[i].erase(model[i].begin());
		NRA = NRA +( getL_NRA(model[i], model[i],G)/ model[i].size());
	}
	y_obj[0] = -NRA;
	//y_graph[1],RC
	double RC = 0;
	for ( int i = 0; i < model.size(); i++){
        if(model[i][0]==model[i][1]&&model[i][1]==0)model[i].erase(model[i].begin());
		RC = RC + (getL_RC(model[i],G) / model[i].size());
	}
	y_obj[1] = RC;
}

void objectives_graph2(vector<double> &x_var, vector <double> &y_obj,MGraph *G){
	vector<vector<double> >model ;
	getModel(x_var,model);
	//y_graph[0],NRA
	double q1 = 0;
	double q2 = 0;
    int model_size=model.size();
    for(int j=0; j < model_size; j++ )
    {
    q1+=(getL_NRA(model[j],model[j],G)/2)/G->numEdges;
    //double q2=(float)(getDegree(model[j],G)*getDegree(model[j],G))/(4*G->numEdges*G->numEdges);
    int kDegree=getDegree(model[j],G);
    q2+=(float)kDegree*kDegree/(4*G->numEdges*G->numEdges);
    }
	y_obj[0] = -q1;
	//y_graph[1],RC
	y_obj[1] = q2;
}
void objectives(vector<double> &x_var, vector <double> &y_obj)
{
	if(!strcmp(strTestInstance,"ZDT1"))
	{
		double g = 0;
		for(int n=1;n<nvar;n++)
			g+= x_var[n];
		g = 1 + 9*g/(nvar-1);

		y_obj[0] = x_var[0];
		y_obj[1] = g*(1 - sqrt(y_obj[0]/g));

        //idle();
	}


	else if(!strcmp(strTestInstance,"ZDT2"))
	{
		double g = 0;
		for(int n=1;n<nvar;n++)
			g+= x_var[n];
		g = 1 + 9*g/(nvar-1);
		y_obj[0] = x_var[0];
		y_obj[1] = g*(1 - pow(y_obj[0]/g,2));

        //idle();
	}

	else if(!strcmp(strTestInstance,"ZDT3"))
	{
		double g = 0;
		for(int n=1;n<nvar;n++)
			g+= x_var[n];
		g = 1 + 9*g/(nvar-1);

		y_obj[0] = x_var[0];
		y_obj[1] = g*(1 - sqrt(x_var[0]/g) - x_var[0]*sin(10*pi*x_var[0])/g);

        //idle();
	}


	else if(!strcmp(strTestInstance,"ZDT4"))
	{
		double g = 0;
		for(int n=1;n<nvar;n++)
		{
			double x = 10*(x_var[n] - 0.5);
			g+= x*x - 10*cos(4*pi*x);
		}
		g = 1 + 10*(nvar-1) + g;
		y_obj[0] = x_var[0];
		y_obj[1] = g*(1- sqrt(y_obj[0]/g));

        //idle();
	}

	else if(!strcmp(strTestInstance,"ZDT6"))
	{
		double g = 0;
		for(int n=1;n<nvar;n++)
			g+= x_var[n]/(nvar - 1);
		g = 1 + 9*pow(g,0.25) ;

		y_obj[0] = 1 - exp(-4*x_var[0])*pow(sin(6*pi*x_var[0]),6);
		y_obj[1] = g*(1- pow(y_obj[0]/g,2));

        //idle();
	}

	// OKA 1
	else if(!strcmp(strTestInstance,"OKA-1"))
	{
		double x1 = 2*pi*(x_var[0] - 0.5);
		double x2 = (x_var[1] - 0.5)*10;
		y_obj[0] = x1;
		y_obj[1] = pi - x1 + fabs(x2 - 5*cos(x1));
	}

	else if(!strcmp(strTestInstance,"OKA-2"))
	{
		double x1 = 2*pow(pi,3)*(x_var[0] - 0.5);
		double x2 = (x_var[1] - 0.5)*10;
		double eta;
		if(x1>=0) eta = pow(x1,1.0/3);
		else      eta = -pow(-x1,1.0/3);

		y_obj[0] = eta;
		y_obj[1] = pi - eta + fabs(x2 - 5*cos(x1));
	}

	else if (!strcmp(strTestInstance, "DTLZ1") )
	{
		DTLZ::DTLZ1(&(*(x_var.begin())),&(*(y_obj.begin())),nvar,nobj);
	}

	else if (!strcmp(strTestInstance, "DTLZ2"))
	{
		DTLZ::DTLZ2(&(*(x_var.begin())), &(*(y_obj.begin())), nvar, nobj);
	}
	else if (!strcmp(strTestInstance, "DTLZ3"))
	{
		DTLZ::DTLZ3(&(*(x_var.begin())), &(*(y_obj.begin())), nvar, nobj);
	}
	else if (!strcmp(strTestInstance, "DTLZ4"))
	{
		DTLZ::DTLZ4(&(*(x_var.begin())), &(*(y_obj.begin())), nvar, nobj);
	}
	else if (!strcmp(strTestInstance, "DTLZ5"))
	{
		DTLZ::DTLZ5(&(*(x_var.begin())), &(*(y_obj.begin())), nvar, nobj);
	}
	else if (!strcmp(strTestInstance, "DTLZ6"))
	{
		DTLZ::DTLZ6(&(*(x_var.begin())), &(*(y_obj.begin())), nvar, nobj);
	}
	else if (!strcmp(strTestInstance, "DTLZ7"))
	{
		DTLZ::DTLZ7(&(*(x_var.begin())), &(*(y_obj.begin())), nvar, nobj);
	}
	else if (!strcmp(strTestInstance, "Convex_DTLZ2"))
	{
		DTLZ::Convex_DTLZ2(&(*(x_var.begin())), &(*(y_obj.begin())), nvar, nobj);
	}
	else if (!strcmp(strTestInstance, "I_DTLZ1"))
	{
		DTLZ::I_DTLZ1(&(*(x_var.begin())), &(*(y_obj.begin())), nvar, nobj);
	}

	else if(!strcmp(strTestInstance,"UF1"))
    {
        CEC09::UF1(&(*(x_var.begin())), &(*(y_obj.begin())),nvar);
    }

    else if(!strcmp(strTestInstance,"UF2"))
    {
        CEC09::UF2(&(*(x_var.begin())), &(*(y_obj.begin())),nvar);
    }
    else if(!strcmp(strTestInstance,"UF3"))
    {
        CEC09::UF3(&(*(x_var.begin())), &(*(y_obj.begin())),nvar);
    }
    else if(!strcmp(strTestInstance,"UF4"))
    {
        CEC09::UF4(&(*(x_var.begin())), &(*(y_obj.begin())),nvar);
    }
    else if(!strcmp(strTestInstance,"UF5"))
    {
        CEC09::UF5(&(*(x_var.begin())), &(*(y_obj.begin())),nvar);
    }
    else if(!strcmp(strTestInstance,"UF6"))
    {
        CEC09::UF6(&(*(x_var.begin())), &(*(y_obj.begin())),nvar);
    }
    else if(!strcmp(strTestInstance,"UF7"))
    {
        CEC09::UF7(&(*(x_var.begin())), &(*(y_obj.begin())),nvar);
    }
    else if(!strcmp(strTestInstance,"UF8"))
    {
        CEC09::UF8(&(*(x_var.begin())), &(*(y_obj.begin())),nvar);
    }
    else if(!strcmp(strTestInstance,"UF9"))
    {
        CEC09::UF9(&(*(x_var.begin())), &(*(y_obj.begin())),nvar);
    }

    else if(!strcmp(strTestInstance,"UF10"))
    {
        CEC09::UF10(&(*(x_var.begin())), &(*(y_obj.begin())),nvar);
    }

}
