#ifndef __COMMON_H_
#define __COMMON_H_


#include "global.h"
//#include <vector>

void minfastsort(vector<double> &x, vector<int> &idx, int n, int m);
void minfastsort(double* x, int* idx, int n, int m);
double dist_vector(vector <double> &vec1, vector <double> &vec2);
double observepoint_dist_vector(vector <double> &vec1, vector <double> &vec2);


template <class T>
void loadpfront(char *filename, vector<T> &ps)
{
    //printf("%s\n",filename);
	std::fstream fin;
	int line=0;
	char str[100]=" ";
	fin.open(filename,std::ios::in);
	if(fin.is_open())
	{
		const char* pch2;
		std::string str;
		while(!fin.eof())
		{
			T  data;
			std::getline(fin,str,'\n');
			pch2 = str.c_str();
			if (pch2 == " ")
			{
				break;
			}
			//==============
			//modified
			//==============
			const char *split = " ";
			char *p = NULL,*pNext = NULL;
			p = strtok_r((char*)pch2, split, &pNext);
			int count_obj = 0;
			while (p != NULL && count_obj < nobj)
			{
				data.y_obj[count_obj] = atof(p);
				count_obj++;
				p = strtok_r(NULL, split, &pNext);
			}

			if (count_obj != nobj)
			{
                //std::cout<<"count_obj:"<<count_obj<<" nobj:"<<nobj<<endl;
				std::cout << "The data of real PF exists errors." << endl;
			}

			line++;
			ps.push_back(data);
		}
	} //end if
	else
		std::cout<<"failed to open "<<filename<<endl;

	fin.close();
}

TCompare ParetoCompare(vector <double> & y_obj1, vector <double>& y_obj2);
void quicksort_increasing(vector <vector<double> > &paretofront, int nPoints);
void save_pop_front(vector <vector<double> >  &pop_front, char saveFilename[1024]);
void random_permutation(int *perm, int size);


//计算组合数
int choose(int n, int m);

double lnchoose(int n,int m);
int cal_index_from_lambda(const vector<int> &V_obj,int n_H);
vector<int> gen_vicinities(vector<int> V_obj,int n_H);

void population2front(vector<vector<double> >&pop_y, vector<vector<double> >&pf_y);
void population2front(vector<vector<double> >&pop_y, vector<vector<double> >&pf_y, vector<int> &pop_index);
void save_population_front(vector <vector<double> >  &pf_y, char saveFilename[1024]);
void save_population_front_split(vector <vector<double> >  &pf_y, vector<int>&pf_index, vector<int>&index_split, char saveFilename[1024]);
void save_pos_split(vector <vector<double> >  &pos, vector<int>&pf_index, vector<int>&index_split, char saveFilename[1024]);
#endif
