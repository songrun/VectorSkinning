#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

typedef Array<double,1,3> Point3d;

int main()
{
	ArrayXXf a(3,3);
	ArrayXXf b(3,3);
	a << 1,2,3,
	   4,5,6,
	   7,8,9;
	b << 1,2,3,
	   1,2,3,
	   1,2,3;
    
    Point3d c;
    Array<double,3,1> d;
    c << 1, 2, 3;
    d << 4, 5, 6;   
       
    a += b; 
	cout << "a + b = " << endl << a << endl;
	
	cout << " c: " << endl << c << endl << c[0] << endl;
	cout << " d: " << endl << d << endl << d[0] << endl;
	
	int u = 5;
	int v = 7;
	float n = floor( min( -log2(u), -log2(v) ) ) + 1;
	cout << "test log " << n << endl;
	
	cout << "test pow " << pow(2, n) << endl;
}