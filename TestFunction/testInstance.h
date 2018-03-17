#ifndef TESTINSTANCE_H
#define TESTINSTANCE_H

#include "../global.h"
#include "cec09.h"
#include "DTLZ.h"

#define pi   3.1415926
#define SQR2  sqrt(2)

void initTestInstance();

void objectives(vector<double> &x_var, vector <double> &y_obj);
double getL_NRA(vector<double> v1, vector<double> v2, MGraph *G);
double getL_RC(vector<double> v1, vector<double> v2, MGraph *G);
vector<vector<double> > getModel(vector<double> x);
void  getModel(vector<double> &x,vector<vector<double> >&model);
void  getModel(vector<double> &x,vector<vector<int> >&model,vector<int>&label);
void getModel(double* x,vector<vector<double> >&model);
void objectives_graph1(vector<double> &x_var, vector <double> &y_obj, MGraph *G);
void objectives_graph2(vector<double> &x_var, vector <double> &y_obj, MGraph *G);
int getDegree(vector<double>& ,MGraph *);

//bool findnum(int l[], int a, int len);
int findmodel(vector<vector<double> >model, int x);
void insertmodel(vector<vector<double> >&model, int p, int index);
void addtemp(vector<double>&temp, int p);

#endif
