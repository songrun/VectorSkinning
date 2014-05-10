/*
clang++ -std=c++11 -stdlib=libc++ -g -I../OpenSubdiv/opensubdiv -I../OpenSubdiv/regression subdivision_matrices_eigen.cpp subdivision_matrices.cpp subdivision_engine.cpp -o subdivision -DDEBUGGING_MAIN -I/opt/local/include/eigen3
*/

#include "subdivision_engine.h"
#include <iostream>
#include <cstdlib>

typedef Eigen::Matrix< subdivision_matrix::real_t, Eigen::Dynamic, 3 > PointList;
typedef Eigen::Matrix< subdivision_matrix::real_t, 1, 3> Point;

const float EPS = 1e-7;
const int resolution = 10;
#define CUBE { {0, 1, 3, 2}, {2, 3, 5, 4}, {4, 5, 7, 6}, {6, 7, 1, 0}, {1, 7, 5, 3}, {6, 0, 2, 4} }


struct prepare_result
{
	Eigen::Matrix< subdivision_matrix::real_t, Eigen::Dynamic, Eigen::Dynamic > inv_system;
	
};

struct subdivision_control_mesh
{
	// This should be an N-by-3 dense matrix.
	PointList vs;
	// This is a vector of faces.
	std::vector< std::vector< int > > faces;
};

Eigen::Matrix< subdivision_matrix::real_t, Eigen::Dynamic, Eigen::Dynamic >
shepard( PointList positions, PointList handle_positions )
{
	using namespace subdivision_matrix;
	
	int row_num = positions.rows();
	int col_num = handle_positions.rows();
	
	Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic > weights( row_num, col_num );
	
	for ( int i = 0; i < row_num; i++ ) {
		int flag = -1;
		for ( int j = 0; j < col_num; j++ ) {
			Point diff = positions.row( i ) - handle_positions.row( j );
			weights( i, j ) = diff.dot( diff );
			// if a position is on the handle, mark the handle's index
			if ( weights( i, j ) <= EPS )
				flag = j;
		}
		// if a position is on any handle, make its weight 1.0.
		if ( flag >= 0 ) {
			for ( int j = 0; j < col_num; j++ ) {
				weights( i, j ) = ( flag == j ? 1.0 : 0.0 );
			}
			continue;
		}
		
		for ( int j = 0; j < col_num; j++ ) {
			weights( i, j ) = 1.0 / weights( i, j );
		}
		weights.row( i ) /= weights.row( i ).sum();
		
	}
	
	return weights;
}


// prepare_result prepare( subdivision_control_mesh, handle_positions, weight_function )
void prepare( const subdivision_control_mesh& mesh, const PointList handle_positions, const std::string weight_function, 
		Eigen::Matrix< subdivision_matrix::real_t, Eigen::Dynamic, Eigen::Dynamic > inv_system )
{
    /// 1 Get the sparse subdivision coefficients at many places around the mesh.
    /// 2 Compute the weights and areas for every point around the mesh.
    
    using namespace subdivision_matrix;
    /// 1
    std::vector< real_t > us, vs;
    createUVs( resolution, resolution, us, vs );
    
    Eigen::SparseMatrix< real_t > M_matrices, Du_matrices, Dv_matrices;
    compute_subdivision_coefficients_for_mesh(
    	mesh.vs.rows(),
     	mesh.faces, 
    	us, vs, 
    	M_matrices, &Du_matrices, &Dv_matrices );
    
    std::cout << M_matrices.rows() << ' ' << M_matrices.cols() << std::endl;
    std::cout << mesh.vs.rows() << ' ' << mesh.vs.cols() << std::endl;
    
    Eigen::Matrix< real_t, Eigen::Dynamic, 3 > positions( M_matrices.rows(), 3 );
    
    positions = M_matrices * mesh.vs;
    
//     std::cout << positions << std::endl;
      
    /// 2
    Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic > weights;
    if( weight_function == "harmonic" ) {
//         weights = harmonic( mesh, positions, handle_positions );

    } else if ( weight_function == "shepard" ) {
        weights = shepard( positions, handle_positions );
    }  
//     std::cout << weights << std::endl;
   
    Eigen::Array< real_t, Eigen::Dynamic, 1 > areas( M_matrices.rows() );
    PointList du_list = Du_matrices * mesh.vs;
    PointList dv_list = Dv_matrices * mesh.vs;
    Eigen::Matrix< real_t, 1, 3 > vec1, vec2;
    for( int i = 0; i < M_matrices.rows(); ++i )
    {
		vec1 = du_list.row( i );
		vec2 = dv_list.row( i );
        areas( i ) = ( vec1.cross( vec2 ).norm() );
    }
//     std::cout << areas << std::endl;
    
//     Eigen::SparseMatrix< real_t > system( M_matrices.cols(), M_matrices.cols() );
	Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic > system( M_matrices.cols(), M_matrices.cols() );
    for( int i = 0; i < 10; ++i )
    {
		Eigen::Matrix< real_t, 1, Eigen::Dynamic > M_i = M_matrices.block( i, 0, 1, M_matrices.cols() );
// 		std::cout << M_i.transpose() * M_i << std::endl; 
   	
        system += M_i.transpose() * M_i * areas( i );
        std::cout << system << std::endl << std::endl;
    }
    
//     return system.lldt(), extra_info;
    return;
}
/*
array solve( prepared, extra_info, handle_transformations )
{
    rhs;
    for( i, handle_transform in handle_transformations )
    {
        rhs += handle_transform * extra_info.positions;
    }
    
    return prepared.solve( rhs );
}*/


void usage( const char* argv0 )
{
    std::cerr << "Usage: " << argv0 << " test_subdivision_engine\n";
    exit( -1 );
}

void print_mesh( const subdivision_control_mesh& mesh );

int main( int argc, char* argv[] )
{
    using std::string;
    
    if( argc != 1 ) usage( argv[0] );
    
	Eigen::Matrix< subdivision_matrix::real_t, Eigen::Dynamic, 3 > handle_positions(2,3);
	handle_positions << 0.8, 0.8, 0.7,
						0.2, -0.3, -0.2;
	
	// make a mesh
	subdivision_control_mesh mesh;
	mesh.vs.resize( 8, 3 );
	mesh.vs <<  0.0,  0.0,  0.0,
				0.0,  0.0,  1.0,
				0.0,  1.0,  0.0,
				0.0,  1.0,  1.0,
				1.0,  0.0,  0.0,
				1.0,  0.0,  1.0,
				1.0,  1.0,  0.0,
				1.0,  1.0,  1.0;
	mesh.faces = CUBE;
	
// 	print_mesh( mesh );
	Eigen::Matrix< subdivision_matrix::real_t, Eigen::Dynamic, Eigen::Dynamic > inv_system;
    prepare( mesh, handle_positions, "shepard", inv_system );

    return 0;
}

void print_mesh( const subdivision_control_mesh& mesh )
{
	std::cout << "vertices\n";
	std::cout << mesh.vs << std::endl;
	std::cout << "faces\n";
	for ( auto face : mesh.faces ) {
		for ( auto v : face )
			std::cout << v << ' ';
		std::cout << std::endl;
	}
}