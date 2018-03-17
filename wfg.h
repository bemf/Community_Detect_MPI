#pragma once
#include "algorithm"
#include <iostream>
#include <vector>
#include <math.h>

using namespace std;

class wfg
{
public:
	wfg(const unsigned int min_max = 0,const unsigned int stop_dimension = 2);
	double compute(vector<vector<double> >& points, const vector<double>& r_point) ;
	vector<double> contributions(std::vector<vector<double> > &, const vector<double> &) ;
	double contributions2(std::vector<vector<double> > &, const vector<double> &);
	~wfg(){}
private:
	void limitset(const unsigned int, const unsigned int, const unsigned int) ;
	double exclusive_hv(const unsigned int, const unsigned int) ;
	double compute_hv(const unsigned int) ;

	///bool cmp_points(double* a, double* b);

	void allocate_wfg_members(std::vector<vector<double> > &, const vector<double> &);
	void free_wfg_members() ;

	const unsigned int m_stop_dimension;
	unsigned int m_max_dim;
	unsigned int m_max_points;
	vector<double> m_refpoint;
	unsigned int m_n_frames;
	unsigned int* m_frames_size;
	double*** m_frames;
};

