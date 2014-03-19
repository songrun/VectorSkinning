#include "SubMatrix.h"
#include <cmath>
#include <Eigen/LU>
#include "cctest_modern.c"


#include <iostream>

SubMatrix::SubMatrix( int n ) : N(n) {

	K = N*2 + 8;
	M = K + 9;
	buildA();
	buildAbar();
	computeEigenvaluesAndEigenvectors();
	computeCoefs();
}

void SubMatrix::buildA() {
	
	A.resize( K, K );
	double a_n = 1-7./(4*N);
	double b_n = 3/( pow(double( N ), 2)*2 );
	double c_n = 1/( pow(double( N ), 2)*4 );
	double a = 9./16;
	double b = 3./32;
	double c = 1./64;
	double d = 3./8;
	double e = 1./16;
	double f = 1./4;
	
	// construct S
	A(0,0) = a_n;
	for ( int i=1; i<2*N+1; i+=2 ) {
		A(0,i) = b_n;
		A(i,0) = d;
		A(0,i+1) = c_n;
		A(i+1,0) = f;
	}
	
	for ( int i=0; i<2*N; i+=2 ) {
		A(1+i, 1+i) = d;
		A(1+i, 1+(i+1)%(2*N)) = e;
		A(1+i, 1+(i+2)%(2*N)) = e;
		A(1+i, 1+(i-1+2*N)%(2*N)) = e;
		A(1+i, 1+(i-2+2*N)%(2*N)) = e;
		
		A(i+2, 2+i) = f;
		A(i+2, 1+i) = f;
		A(i+2, 1+(i+2)%(2*N)) = f;
	}	
	
	// construct S11
	if( N == 3 ) {
		A.bottomLeftCorner(7, 7) << 
			c, 0, 0, b, a, b, 0,
			e, 0, 0, e, d, d, 0,
			b, c, 0, c, b, a, b,
			e, e, 0, 0, 0, d, d,
			e, 0, 0, d, d, e, 0,
			b, c, b, a, b, c, 0,
			e, e, d, d, 0, 0, 0;
	}
	else {
		A.bottomLeftCorner(7, 9) << 
			c, 0, 0, b, a, b, 0, 0, 0,
			e, 0, 0, e, d, d, 0, 0, 0,
			b, 0, 0, c, b, a, b, c, 0,
			e, 0, 0, 0, 0, d, d, e, 0,
			e, 0, 0, d, d, e, 0, 0, 0,
			b, c, b, a, b, c, 0, 0, 0,
			e, e, d, d, 0, 0, 0, 0, 0;
	}
	
	// construct S12
	A.bottomRightCorner(7, 7) <<
		c, b, c, 0, b, c, 0,
		0, e, e, 0, 0, 0, 0,
		0, c, b, c, 0, 0, 0,
		0, 0, e, e, 0, 0, 0,
		0, 0, 0, 0, e, e, 0,
		0, 0, 0, 0, c, b, c,
		0, 0, 0, 0, 0, e, e;
			
}

void SubMatrix::buildAbar() {

	Abar.resize( M, K );
 	Abar.block(0, 0, K, K) = A.block(0, 0, K, K);
	
	double d = 3./8;
	double e = 1./16;
	double f = 1./4;
	
	Abar.bottomLeftCorner(9, 7) << 
		0, 0, 0, 0, f, 0, 0,
		0, 0, 0, 0, d, e, 0,
		0, 0, 0, 0, f, f, 0,
		0, 0, 0, 0, e, d, e,
		0, 0, 0, 0, 0, f, f,
		0, 0, 0, e, d, 0, 0,
		0, 0, 0, f, f, 0, 0,
		0, 0, e, d, e, 0, 0,
		0, 0, f, f, 0, 0, 0;
		
	Abar.bottomRightCorner(9, 7) <<
		f, f, 0, 0, f, 0, 0,
		e, d, e, 0, e, 0, 0,
		0, f, f, 0, 0, 0, 0,
		0, e, d, e, 0, 0, 0,
		0, 0, f, f, 0, 0, 0,
		e, e, 0, 0, d, e, 0,
		0, 0, 0, 0, f, f, 0,
		0, 0, 0, 0, e, d, e,
		0, 0, 0, 0, 0, f, f;		
}

void SubMatrix::computeEigenvaluesAndEigenvectors() {
	es.compute( A, true );
}

void SubMatrix::computeCoefs() {

	MatrixXd Vbar = Abar * (es.eigenvectors().real());
	coefs.resize( K, 16 );
	Matrix4d BN;
	BN <<  	1, -3, 3, -1,
			4, 0, -6, 3,
			1, 3, 3, -3, 
			0, 0, 0, 1;
	BN /= 6.0;
	
	VectorXi q1(16);
	Matrix4d temp;
	q1 << 8, 7, 2*N+5, 2*N+13, 1, 6, 2*N+4, 2*N+12, 4, 5, 2*N+3, 2*N+11, 2*N+7, 2*N+6, 2*N+2, 2*N+10; 
	for( int i=0; i<K; i++ ) {
		for ( int j=0; j<16; j++ ) {

 			temp = ( BN.row(j%4).transpose() * BN.row(j/4) ) * Vbar(q1(j), i);
 			for ( int t=0; t<4; t+=4 ) {
  				coefs.block<1,4>(i,t) += temp.row(t);
  			}
		}
	}
	
}

void SubMatrix::print( ) {
	cout << " matrix size: " << A.size() << endl;
	cout << A << endl;
	
	cout << " Abar matrix size: " << Abar.size() << endl;
	cout << Abar << endl;
	
	cout << "The eigenvalues of A are:" << endl << es.eigenvalues().real() << endl;
	cout << "The matrix of eigenvectors, V, is:" << endl << es.eigenvectors().real() << endl << endl;
	
	cout << "Vbar equals to Abar * V, which is:" << endl << Abar * (es.eigenvectors().real()) << endl << endl;
	
	cout << "Coefficients matrix is:" << endl << coefs << endl << endl;
// 	for (int i=0; i<A.cols(); i++) {
// 		cout << "A*v1: " << A * es.eigenvectors().real().col(i) << endl;
// 		cout << "alpha1*v1: " << es.eigenvalues().real()[i]* es.eigenvectors().real().col(i) << endl << endl;
// 	}

// 	MatrixXd ref(14, 14);
// 	
// 	ref << 
// 	 0.375,  0.167,  0.042,  0.167,  0.042,  0.167,  0.042,  0.000, -0.000,  0.000, -0.000, -0.000,  0.000, -0.000, 
// 	-0.000,  0.000, 0.219,  0.560, -0.000, -0.560, -0.219,  0.000, -0.000,  0.000, -0.000, -0.000,  0.000,  0.000, 
// 	 0.000,  0.647,  0.126, -0.323, -0.252, -0.323,  0.126,  0.000, -0.000,  0.000,  0.000, -0.000,  0.000, -0.000, 
// 	-0.450, -0.000,  0.150, -0.000,  0.150, -0.000,  0.150, -0.000,  0.000,  0.000, -0.000,  0.000, -0.000,  0.000, 
// 	 0.000, -0.000,  0.359, -0.560, -0.000,  0.560, -0.359, -0.000,  0.000, -0.000,  0.000,  0.000,  0.000, -0.000, 
// 	-0.000, -0.647,  0.207,  0.323, -0.414,  0.323,  0.207, -0.000,  0.000, -0.000,  0.000, -0.000,  0.000, -0.000, 
// 	 3.000,  1.333, -1.167, -2.333, -1.167,  1.333, -2.000, -0.000,  0.000, -0.000,  0.000,  0.167,  0.667,  0.167, 
// 	 3.000,  1.333, -2.000,  1.333, -1.167, -2.333, -1.167,  0.000,  0.167,  0.667,  0.167, -0.000,  0.000, -0.000, 
// 	-0.000, -1.312,  0.000,  1.312, -1.437,  0.000,  1.437,  0.000,  0.500,  0.000, -0.500,  0.000, -0.000,  0.000, 
// 	 0.000, -1.312,  1.437, -0.000, -1.437, 1.312,  0.000, -0.000, -0.000, -0.000,  0.000,  0.500, -0.000, -0.500, 
// 	 1.615, -0.720, -0.077, -0.227, -0.077, -0.720,  0.206, -0.000,  0.000, -0.000,  0.000,  0.167, -0.333,  0.167, 
// 	 1.615, -0.720,  0.206, -0.720, -0.077, -0.227, -0.077, -0.000,  0.167, -0.333,  0.167,  0.000, -0.000,  0.000, 
// 	-41.897, 19.683, -3.659, 14.158,  2.214, 14.158, -3.659,  1.000, -3.000,  3.000, -1.000, -3.000,  3.000, -1.000, 
// 	 0.375, -0.167,  0.042, -0.167,  0.042, -0.167,  0.042, -0.000,  0.000, -0.000,  0.000,  0.000, -0.000,  0.000;
// 	cout << "The invserse of The matrix of eigenvectors, V, is:" << endl << ref.inverse() << endl << endl;

}

int main( int argc, char ** argv ) {
	
	SubMatrix A = SubMatrix(3);
	A.print();

}
