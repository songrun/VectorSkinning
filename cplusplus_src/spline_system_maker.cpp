/*
	A C++ version of the python code generate_chain_system.py.
*/

#include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

struct Bundle {
	
	Eigen::Matrix4f [] W_matrices;
	Eigen::Matrix43f control_points;
	Eigen::Array22f constraints;
	float length;
	Eigen::ArrayXf ts;
	Eigen::ArrayXf dts;
	Eigen::Array2f magnitudes;
	Eigen::Array22f directions; 
	
	Bundle( Eigen::Matrix4f [] a, Eigen::Matrix43f b, Eigen::Array22f c, float d, Eigen::ArrayXf e, Eigen::ArrayXf f, Eigen::Array2f [] g, Eigen::Array22f [] h ) : 
	W_matrices(a), control_points(b), constraints(c), length(d), ts(e), dts(f), magnitudes(g), directions(h) 
	{
		if magnitudes == NULL:
		for( int i=0; i<control_points.size(); i++ ) {
			magnitudes[0] = control_points.row(i) - control_points.row(i)
		}	
	}
}



class BezierConstraintSolver {

public:
	void build_system
	
private:
	Bundle [] bundles;
}

BezierConstraintSolver::BezierConstraintSolver () : 