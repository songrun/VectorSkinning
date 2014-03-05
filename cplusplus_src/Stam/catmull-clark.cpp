#include <iostream>
#include <Eigen/Dense>
#include "SubMatrix.h"
#include "cctest_modern.c"

using namespace std;
using namespace Eigen;

typedef Array<double,1,3> Point3d;

typedef
   struct
   {
      VectorXd values;
      MatrixXd vectorI;
      Array<MatrixXd,3,1>  coefs;
   } EIGENSPACE;
   
EIGENSPACE * MakeEigenSpace( EVALSTRUCT * ev ) 
{
	
}

/*
///  Step 1:
void ProjectPoints( EVALSTRUCT ** ev, Point3d *Cp, Point3d *C, int N ) {
	
	for ( int i=0; i<2*N+8; i++ ) {
		Cp[i] << 0, 0, 0;
		for ( int j=0; j<2*N+8; j++ ) {
			Cp[i] += ev[ N ].iV[i][j] * C[j];
		}
	}
}

///  Step 2:
void EvalSurf( EVALSTRUCT ** ev, Point3d P, double u, double v, Point3d *Cp, int N ) {
	// determine in which domain Omega(n, k) the parameter lies
	n = floor( min( -log2(u), -log2(v) ) ) + 1;
	pow2 = pow(2, n-1);
	u *= pow2;
	v *= pow2;
	
	if ( v < 0.5 ) {
		k = 0;
		u = 2*u-1;
		v = 2*v;
	}
	else if ( u < 0.5 ) {
		k = 2;
		u = 2*u;
		v = 2*v-1;
	}
	eles {
		k = 1;
		u = 2*u-1;
		v = 2*v-1;
	}
	
	// Evaluate the surface
	P << 0, 0, 0;
	for ( int i=0; i<2*N+8; i++ ) {
		P += pow( eigen[N].val[i], n-1 ) * EvalSpline( ev[N].Phi[i][k], u, v ) * Cp[i];
	}
}

double EvalUniformCubicSplineBasis( int i, double t ) {
	
	double basis = -1;
	
	switch( i ) {
		case 0:
			basis = ( 1 - 3*t + 3*t*t - t*t*t )/6;
		case 1:
			basis = 4 - 6*t*t + 3*t*t*t;
		case 2:
			basis = 1 + 3*t + 3*t*t -3*t*t*t;
		case 3:
			basis = t*t*t;
		default:	
			cerr << "unexpected index for cubic spline basis.";
	}
	
	return basis;
}

double EvalSpline( double *coef, double u, double v ) {
	
	Point3d p;
	
	for ( int i=0; i<16; i++ ) {
		
	}
} 
*/

int main() {
	
	int Nmax;
	EVALSTRUCT ** ev;
	
   	ev = read_eval ( &Nmax );
   	if ( !ev ) exit ( 1 );
	
	print_eval_N ( ev, 3 );
	
	EIGENSPACE * es = MakeEigenSpace( ev[3] );
	
	return 0;
}