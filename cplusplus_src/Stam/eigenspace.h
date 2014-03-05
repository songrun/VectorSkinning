#include <iostream>
#include "definition.h"
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

class EigenSpace {

public:
	EigenSpace( int );
	void buildA();
	void buildAbar();
	void computeCoefs();
	void print();

private:
	MatrixXd A;
	MatrixXd Abar;
	EigenSolver<MatrixXd> es;
	MatrixXd coefs;
	int N;
	int K;
	int M;
};