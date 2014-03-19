#include <iostream>
#include <Eigen/Dense>
#include "SubMatrix.h"
#include "cctest_modern.c"

using namespace std;

typedef Array<double,1,3> Point3d;

typedef
   struct
   {
      VectorXd values;
      MatrixXd vectorI;
      Array<MatrixXd,3,1>  coefs;
   } EIGENSPACE;
   
EIGENSPACE * makeEigenSpace( EVALSTRUCT ** ev, int Nmax ) 
{
	EIGENSPACE * es;
	es = new EIGENSPACE [Nmax - 2];
	
	for ( int p=0; p<=Nmax-3; p++ )
	{
		int N = p+3;
		int K = 2*N+8;
	
		es -> values = VectorXd( K );
		for ( int i=0; i<K; i++ ) {
			es -> values(i) = ev[p] -> val[i];
		}
	
		es -> vectorI = MatrixXd( K, K );
		for ( int i=0; i<K; i++ ) {
			for ( int j=0; j<K; j++ ) {
				es -> vectorI(i, j) = ev[p] -> vecI[i+K*j];
			}
		}
	
		for ( int k=0; k<3; k++ ) {
			es -> coefs(k) = MatrixXd( K, 16 ); 
			for ( int i=0; i<K; i++ ) {
				for ( int j=0; j<16; j++ ) {
					es -> coefs(k)(i,j) = ev[p] -> Phi[k][i+K*j];
				}
			}
		}
	}
	return es;
	
}

void printEigenSpace( EIGENSPACE *p, int N ) 
{
	EIGENSPACE *es = p+N-2;
	
	cout << "Eigen eigenvalues: " << endl << es -> values << endl;
	cout << "Eigen inverse of eigenvecters: " << endl << es -> vectorI << endl;
	for ( int k=0; k<3; k++)
 	cout << "Eigen coefficients: " << k << endl << es -> coefs(k) << endl;
}

/*
///  Step 1:
void projectPoints( EIGENSPACE ** ev, Point3d *Cp, Point3d *C, int N ) {
	
	for ( int i=0; i<2*N+8; i++ ) {
		Cp[i] << 0, 0, 0;
		for ( int j=0; j<2*N+8; j++ ) {
			Cp[i] += ev[ N ].iV[i][j] * C[j];
		}
	}
}
/*
///  Step 2:
void evalSurf( EVALSTRUCT ** ev, Point3d P, double u, double v, Point3d *Cp, int N ) {
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

double evalUniformCubicSplineBasis( int i, double t ) {
	
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

double evalSpline( double *coef, double u, double v ) {
	
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
	
// 	print_eval_N ( ev, 3 );
	
	EIGENSPACE * es = makeEigenSpace( ev, 3 );
 	printEigenSpace( es, 3 );
	
	return 0;
}