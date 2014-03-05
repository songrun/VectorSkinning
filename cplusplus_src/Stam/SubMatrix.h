#include <iostream>
#include "definition.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

using namespace std;
using namespace Eigen;

class SubMatrix {

public:
	SubMatrix( int );
	void buildA();
	void buildAbar();
	void computeEigenvaluesAndEigenvectors();
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
