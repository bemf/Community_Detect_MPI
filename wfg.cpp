#include "wfg.h"

unsigned int m_current_slice;
unsigned int min_max;//min 0,max = 1;

// Domination results of the 'dom_cmp' methods
enum {
	DOM_CMP_B_DOMINATES_A = 1, // second argument dominates the first one
	DOM_CMP_A_DOMINATES_B = 2, // first argument dominates the second one
	DOM_CMP_A_B_EQUAL = 3, // both points are equal
	DOM_CMP_INCOMPARABLE = 4 // points are incomparable
};

struct Less{
    bool operator()(double* a, double* b)
    {
        for (int i = m_current_slice - 1; i >= 0; --i){
            if (a[i] > b[i]) {
                return true;
            }
            else if (a[i] < b[i]) {
                return false;
            }
        }
        return false;

    }
};

struct LessOne{
    bool operator()(double* a, double* b)
    {
      //  for (int i = m_current_slice - 1; i >= 0; --i)
		if(m_current_slice - 1 >=0){
            if (a[ m_current_slice - 1] > b[ m_current_slice - 1]) {
                return true;
            }
            else if (a[ m_current_slice - 1] < b[ m_current_slice - 1]) {
                return false;
            }
        }
        return false;

    }
};

/// Dominance comparison method
/**
* Establishes the domination relationship between two points.
*
* returns DOM_CMP_B_DOMINATES_A if point 'b' DOMINATES point 'a'
* returns DOM_CMP_A_DOMINATES_B if point 'a' DOMINATES point 'b'
* returns DOM_CMP_A_B_EQUAL if point 'a' IS EQUAL TO 'b'
* returns DOM_CMP_INCOMPARABLE otherwise
*/
int dom_cmp(double* &a, double* &b, unsigned int dim_bound)
{
	if (dim_bound == 0) {
		dim_bound = sizeof(a) / sizeof(double);
	}
	if(!min_max){///最小化问题
		for (int i = 0; i < dim_bound; ++i) {
			if (a[i] > b[i]) {
				for (int j = i + 1; j < dim_bound; ++j) {
					if (a[j] < b[j]) {
						return DOM_CMP_INCOMPARABLE;
					}
				}
				return DOM_CMP_B_DOMINATES_A;
			}
			else if (a[i] < b[i]) {
				for (int j = i + 1; j < dim_bound; ++j) {
					if (a[j] > b[j]) {
						return DOM_CMP_INCOMPARABLE;
					}
				}
				return DOM_CMP_A_DOMINATES_B;
			}
		}
		return DOM_CMP_A_B_EQUAL;
	}
	else{///最大化问题
		for (int i = 0; i < dim_bound; ++i) {
			if (a[i] < b[i]) {
				for (int j = i + 1; j < dim_bound; ++j) {
					if (a[j] > b[j]) {
						return DOM_CMP_INCOMPARABLE;
					}
				}
				return DOM_CMP_B_DOMINATES_A;
			}
			else if (a[i] > b[i]) {
				for (int j = i + 1; j < dim_bound; ++j) {
					if (a[j] < b[j]) {
						return DOM_CMP_INCOMPARABLE;
					}
				}
				return DOM_CMP_A_DOMINATES_B;
			}
		}
		return DOM_CMP_A_B_EQUAL;
	}
}
/// Compute volume between two points
/**
* Calculates the volume between points a and b (as defined for n-dimensional Euclidean spaces).
*
* @param[in] a first point defining the hypercube
* @param[in] b second point defining the hypercube
* @param[in] dim_bound dimension boundary for the volume. If equal to 0 (default value), then compute the volume of whole vector. Any positive number limits the computation from dimension 1 to dim_bound INCLUSIVE.
*
* @return volume of hypercube defined by points a and b
*/
double volume_between(double* &a, vector<double> &b, unsigned int dim_bound = 0)
{
	if (dim_bound == 0) {
		dim_bound = sizeof(a) / sizeof(double);
	}
	double volume = 1.0;
	for (int idx = 0; idx < dim_bound; ++idx) {
		volume *= ( b[idx] - a[idx]);
	}
	return fabs(volume);//(volume < 0 ? -volume : volume);
}


/// Constructor
wfg::wfg(const unsigned int c_min_max ,const unsigned int stop_dimension):m_stop_dimension(stop_dimension)
{
    m_current_slice = 0;
	min_max = c_min_max;
	if (stop_dimension < 2) {
		cout << "Stop dimension for WFG must be greater than or equal to 2" << endl;
	}
}

/// Compute hypervolume
/**
* Computes the hypervolume using the WFG algorithm.
*
* @param[in] points vector of points containing the D-dimensional points for which we compute the hypervolume
* @param[in] r_point reference point for the points
*
* @return hypervolume.
*/
double wfg::compute(std::vector<vector<double> > &points, const vector<double> &r_point)
{
	vector<vector<double> > filterPoints;
	//pre work
	for(int i=0;i<points.size();i++)
	{
		bool isAdd = true;
		for(int j=0;j<points[i].size();j++){
			if(!min_max && points[i][j] > r_point[j])
			{
				isAdd = false;
				break;
			}
			else if(min_max && points[i][j] < r_point[j])
			{
				isAdd = false;
				break;
			}
		}
		if(isAdd)
			filterPoints.push_back(points[i]);
	}
	if(filterPoints.size() == 0)
		return 0;
	
	allocate_wfg_members(filterPoints, r_point);
	double hv = compute_hv(1);
	free_wfg_members();
	return hv;
}

/// Contributions method
/**
* This method employs a slightly modified version of the original WFG algorithm to suit the computation of the exclusive contributions.
* It differs from the IWFG algorithm (referenced below), as we do not use the priority-queueing mechanism, but compute every exclusive contribution instead.
* This may suggest that the algorithm for the extreme contributor itself reduces to the 'naive' approach. It is not the case however,
* as we utilize the benefits of the 'limitset', before we begin the recursion.
* This simplifies the sub problems for each exclusive computation right away, which makes the whole algorithm much faster, and in many cases only slower than regular WFG algorithm by a constant factor.
*
* @see "Lyndon While and Lucas Bradstreet. Applying the WFG Algorithm To Calculate Incremental Hypervolumes. 2012 IEEE Congress on Evolutionary Computation. CEC 2012, pages 489-496. IEEE, June 2012."
*
* @param[in] points vector of points containing the D-dimensional points for which we compute the hypervolume
* @param[in] r_point reference point for the points
*/
vector<double> wfg::contributions(std::vector<vector<double> > &points, const vector<double> &r_point)
{


	std::vector<double> c;
	c.reserve(points.size());

	// Allocate the same members as for 'compute' method
	allocate_wfg_members(points, r_point);

	// Prepare the memory for first front
	double** fr = new double*[m_max_points];
	for (unsigned int i = 0; i < m_max_points; ++i) {
		fr[i] = new double[m_current_slice];
	}
	m_frames[m_n_frames] = fr;
	m_frames_size[m_n_frames] = 0;
	++m_n_frames;

	for (unsigned int p_idx = 0; p_idx < m_max_points; ++p_idx) {
		limitset(0, p_idx, 1);
		c.push_back(exclusive_hv(p_idx, 1));
	}

	// Free the contributions and the remaining WFG members
	free_wfg_members();

	return c;
}


double wfg::contributions2(std::vector<vector<double> > &points, const vector<double> &r_point)
{
	
	// Allocate the same members as for 'compute' method
	allocate_wfg_members(points, r_point);

	// Prepare the memory for first front
	double** fr = new double*[m_max_points];
	for (unsigned int i = 0; i < m_max_points; ++i) {
		fr[i] = new double[m_current_slice];
	}
	m_frames[m_n_frames] = fr;
	m_frames_size[m_n_frames] = 0;
	++m_n_frames;

	limitset(0, 0, 1);
	double result = exclusive_hv(0, 1);
	// Free the contributions and the remaining WFG members
	free_wfg_members();

	return result;
}

/// Allocate the memory for the 'compute' method
void wfg::allocate_wfg_members(std::vector<vector<double> > &points, const vector<double> &r_point)
{
	m_max_points = points.size();
	m_max_dim = r_point.size();

	m_refpoint = r_point;

	// Reserve the space beforehand for each level or recursion.
	// WFG with slicing feature will not go recursively deeper than the dimension size.
	m_frames = new double**[m_max_dim];
	m_frames_size = new unsigned int[m_max_dim];

	// Copy the initial set into the frame at index 0.
	double** fr = new double*[m_max_points];
	for (unsigned int p_idx = 0; p_idx < m_max_points; ++p_idx) {
		fr[p_idx] = new double[m_max_dim];
		for (unsigned int d_idx = 0; d_idx < m_max_dim; ++d_idx) {
			fr[p_idx][d_idx] = points[p_idx][d_idx];
		}
	}
	m_frames[0] = fr;
	m_frames_size[0] = m_max_points;
	m_n_frames = 1;

	// Variable holding the current "depth" of dimension slicing. We progress by slicing dimensions from the end.
	m_current_slice = m_max_dim;
}

/// Free the previously allocated memory
void wfg::free_wfg_members()
{
	// Free the memory.
	m_refpoint.clear();

	for (unsigned int fr_idx = 0; fr_idx < m_n_frames; ++fr_idx) {
		for (unsigned int p_idx = 0; p_idx < m_max_points; ++p_idx) {
			delete[] m_frames[fr_idx][p_idx];
		}
		delete[] m_frames[fr_idx];
	}
	delete[] m_frames;
	delete[] m_frames_size;
}

/// Limit the set of points to point at p_idx
void wfg::limitset(const unsigned int begin_idx, const unsigned int p_idx, const unsigned int rec_level)
{
	double **points = m_frames[rec_level - 1];
	unsigned int n_points = m_frames_size[rec_level - 1];

	int no_points = 0;

	double* p = points[p_idx];
	double** frame = m_frames[rec_level];

	for (unsigned int idx = begin_idx; idx < n_points; ++idx) {
		if (idx == p_idx) {
			continue;
		}
		
		if(min_max){//最大化
			for (int f_idx = 0; f_idx < m_current_slice; ++f_idx) {
				frame[no_points][f_idx] = std::min(points[idx][f_idx], p[f_idx]);
			}	
		}
		else{//最小化
			for (int f_idx = 0; f_idx < m_current_slice; ++f_idx) {
				frame[no_points][f_idx] = std::max(points[idx][f_idx], p[f_idx]);
			}
		}

		std::vector<int> cmp_results;
		cmp_results.resize(no_points);
		double* s = frame[no_points];

		bool keep_s = true;

		// Check whether any point is dominating the point 's'.
		for (int q_idx = 0; q_idx < no_points; ++q_idx) {
			cmp_results[q_idx] = dom_cmp(s, frame[q_idx], m_current_slice);
			if (cmp_results[q_idx] == DOM_CMP_B_DOMINATES_A) {
				keep_s = false;
				break;
			}
		}
		// If neither is, remove points dominated by 's' (we store that during the first loop).
		if (keep_s) {
			int prev = 0;
			int next = 0;
			while (next < no_points) {
				if (cmp_results[next] != DOM_CMP_A_DOMINATES_B && cmp_results[next] !=DOM_CMP_A_B_EQUAL) {
					if (prev < next) {
						for (unsigned int d_idx = 0; d_idx < m_current_slice; ++d_idx) {
							frame[prev][d_idx] = frame[next][d_idx];
						}
					}
					++prev;
				}
				++next;
			}
			// Append 's' at the end, if prev==next it's not necessary as it's already there.
			if (prev < next) {
				for (unsigned int d_idx = 0; d_idx < m_current_slice; ++d_idx) {
					frame[prev][d_idx] = s[d_idx];
				}
			}
			no_points = prev + 1;
		}
	}

	m_frames_size[rec_level] = no_points;
}

/// Compute the hypervolume recursively
double wfg::compute_hv(const unsigned int rec_level)
{
	
	double **points = m_frames[rec_level - 1];
	unsigned int n_points = m_frames_size[rec_level - 1];

	// Simple inclusion-exclusion for one and two points
	if (n_points == 1) {
		return volume_between(points[0], m_refpoint, m_current_slice);;
	}
	else if (n_points == 2) {
		
		double hv = volume_between(points[0], m_refpoint, m_current_slice)
			+ volume_between(points[1], m_refpoint, m_current_slice);
		double isect = 1.0;
		for (unsigned int i = 0; i<m_current_slice; ++i) {
			if(!min_max)///最小化
				isect *= (m_refpoint[i] - std::max(points[0][i], points[1][i]));
			else//最大化
				isect *= (std::min(points[0][i],points[1][i]) - m_refpoint[i]);
		}

		return hv - isect;//(isect < 0 ? -isect : isect);
	}

	// If already sliced to dimension at which we use another algorithm.
	if (m_current_slice == m_stop_dimension) {
		
		if (m_stop_dimension == 2) {
			// Use a very efficient version of hv2d
			if (n_points == 0) {
				return 0.0;
			}
			else if (n_points == 1) {
				return volume_between(points[0], m_refpoint);
			}

			sort(points, points + n_points,Less());

			double hypervolume = 0.0;

			// width of the sweeping line
			double w = m_refpoint[0] - points[0][0];
			if(min_max){//最大化
				for (unsigned int idx = 0; idx < n_points - 1; ++idx) {
					hypervolume += (points[idx][1] - points[idx + 1][1]) * w;
					w = std::max(w,points[idx + 1][0] - m_refpoint[0]);
				}
				hypervolume += (points[n_points - 1][1] - m_refpoint[1]) * w;
			}
			else{//最小化
				hypervolume += (m_refpoint[1] - points[0][1]) * w;
				for (unsigned int idx = 0; idx < n_points - 1; ++idx) {
					hypervolume += (points[idx][1] - points[idx + 1][1]) * w;
					w = std::min(w, m_refpoint[0] - points[idx + 1][0]);
				}
			}
			return hypervolume;
		}
		else {
			// Let hypervolume object pick the best method otherwise.
			/*
			std::vector<fitness_vector> points_cpy;
			points_cpy.reserve(n_points);
			for (unsigned int i = 0; i < n_points; ++i) {
				points_cpy.push_back(fitness_vector(points[i], points[i] + m_current_slice));
			}
			fitness_vector r_cpy(m_refpoint, m_refpoint + m_current_slice);

			hypervolume hv = hypervolume(points_cpy, false);
			hv.set_copy_points(false);
			return hv.compute(r_cpy);
			*/
		}
	}
	else {
		// Otherwise, sort the points in preparation for the next recursive step
		// Bind the object under "this" pointer to the cmp_points method so it can be used as a valid comparator function for std::sort
		// We need that in order for the cmp_points to have acces to the m_current_slice member variable.
		sort(points, points + n_points,LessOne());
	}

	double H = 0.0;
	--m_current_slice;

	if (rec_level >= m_n_frames) {
		double** fr = new double*[m_max_points];
		for (unsigned int i = 0; i < m_max_points; ++i) {
			fr[i] = new double[m_current_slice];
		}
		m_frames[m_n_frames] = fr;
		m_frames_size[m_n_frames] = 0;
		++m_n_frames;
	}

	for (unsigned int p_idx = 0; p_idx < n_points; ++p_idx) {
		
		limitset(p_idx + 1, p_idx, rec_level);

		H += (fabs(m_refpoint[m_current_slice] - points[p_idx][m_current_slice]) * exclusive_hv(p_idx, rec_level));
	
	}
	++m_current_slice;
	return H;
}

/// Compute the exclusive hypervolume of point at p_idx
double wfg::exclusive_hv(const unsigned int p_idx, const unsigned int rec_level)
{
	//double H = base::volume_between(points[p_idx], m_refpoint, m_current_slice);
	double H =volume_between(m_frames[rec_level - 1][p_idx], m_refpoint, m_current_slice);

	if (m_frames_size[rec_level] == 1) {
		H -= volume_between(m_frames[rec_level][0], m_refpoint, m_current_slice);
	}
	else if (m_frames_size[rec_level] > 1) {
		H -= compute_hv(rec_level + 1);
	}

	return H;
}


