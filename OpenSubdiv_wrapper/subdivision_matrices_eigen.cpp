// c++ -g -std=c++11 -I/Users/yotam/Work/ext/OpenSubdiv/opensubdiv -I/Users/yotam/Work/ext/OpenSubdiv/regression subdivision_matrices_eigen.cpp subdivision_matrices.cpp -o subdivision_matrices -DDEBUGGING_MAIN -I/usr/local/include/eigen3
// clang++ -std=c++11 -g -I/Users/Solomon/Graphics/OpenSubdiv/opensubdiv -I/Users/Solomon/Graphics/OpenSubdiv/regression subdivision_matrices_eigen.cpp subdivision_matrices.cpp -o subdivision_matrices -DDEBUGGING_MAIN -I/opt/local/include/eigen3

#include "subdivision_matrices_eigen.h"

#include <cassert>

namespace
{
/*
Given the number of vertices 'num_vertices',
a vector of sparse_vector_t 'sparse_vectors' containing a sequence of ( index, coefficient ) pairs,
and an output Eigen::SparseMatrix 'matrix',
fills 'matrix' such that the matrix multiplication of 'matrix' times a num_vertices-by-K matrix of control points
yields the positions specified by 'sparse_vectors'.
*/
template< typename T >
void convert_vector_of_sparse_vectors_to_matrix( int num_vertices, const std::vector< subdivision_matrix::sparse_vector_t >& sparse_vectors, Eigen::SparseMatrix< T >& matrix )
{
    assert( num_vertices > 0 );
    
    // Clear the output matrix.
    matrix.resize( sparse_vectors.size(), num_vertices );
    matrix.setZero();
    
    // We will fill the matrix with a vector of triplets.
    std::vector< Eigen::Triplet< T > > triplets;
    
    // Count the number of triplets we will need.
    int nnz = 0;
    for( const auto sp : sparse_vectors ) nnz += sp.size();
    triplets.reserve( nnz );
    
    // Convert 'sparse_vectors' to triplets.
    for( unsigned int row = 0; row < sparse_vectors.size(); ++row )
    {
        for( const auto p : sparse_vectors[row] )
        {
            triplets.push_back( { row, p.first, p.second } );
        }
    }
    
    matrix.setFromTriplets( triplets.begin(), triplets.end() );
}
}

namespace subdivision_matrix
{

void compute_subdivision_coefficients_for_mesh(
    int num_vertices,
    const std::vector< std::vector< int > >& faces,
    const std::vector< real_t >& us,
    const std::vector< real_t >& vs,
    Eigen::SparseMatrix< real_t >& positions_out,
    Eigen::SparseMatrix< real_t >* du_out,
    Eigen::SparseMatrix< real_t >* dv_out
    )
{
    // Intermediary output vectors.
    std::vector< sparse_vector_t > position_coeffs;
    std::vector< sparse_vector_t > du_coeffs;
    std::vector< sparse_vector_t > dv_coeffs;
    
    // Call the function we're wrapping.
    compute_subdivision_coefficients_for_mesh(
        num_vertices,
        faces,
        us, vs,
        position_coeffs,
        du_out ? &du_coeffs : nullptr,
        dv_out ? &dv_coeffs : nullptr
        );
    
    convert_vector_of_sparse_vectors_to_matrix( num_vertices, position_coeffs, positions_out );
    if( du_out )
    {
        assert( dv_out );
        
        convert_vector_of_sparse_vectors_to_matrix( num_vertices, du_coeffs, *du_out );
        convert_vector_of_sparse_vectors_to_matrix( num_vertices, dv_coeffs, *dv_out );
    }
}

}

// test
void print_sparse_vectors( const std::vector< subdivision_matrix::sparse_vector_t >& sparse_vectors, const std::string& label );

#define TORUS 32, { {4, 5, 1, 0}, {5, 6, 2, 1}, {6, 7, 3, 2}, {7, 4, 0, 3}, {8, 9, 5, 4}, {9, 10, 6, 5}, {10, 11, 7, 6}, {11, 8, 4, 7}, {12, 13, 9, 8}, {13, 14, 10, 9}, {14, 15, 11, 10}, {15, 12, 8, 11}, {16, 17, 13, 12}, {17, 18, 14, 13}, {18, 19, 15, 14}, {19, 16, 12, 15}, {20, 21, 17, 16}, {21, 22, 18, 17}, {22, 23, 19, 18}, {23, 20, 16, 19}, {24, 25, 21, 20}, {25, 26, 22, 21}, {26, 27, 23, 22}, {27, 24, 20, 23}, {28, 29, 25, 24}, {29, 30, 26, 25}, {30, 31, 27, 26}, {31, 28, 24, 27}, {0, 1, 29, 28}, {1, 2, 30, 29}, {2, 3, 31, 30}, {3, 0, 28, 31} }
#define CUBE 8, { {0, 1, 3, 2}, {2, 3, 5, 4}, {4, 5, 7, 6}, {6, 7, 1, 0}, {1, 7, 5, 3}, {6, 0, 2, 4} }


#include <iostream>
#include <cstdlib>
void test_compute_coefficients( int resolution, Eigen::SparseMatrix< subdivision_matrix::real_t >& position_out )
{
    using namespace subdivision_matrix;
    
    std::vector< real_t > us, vs;
    createUVs( resolution, resolution, us, vs );
    
//   	Eigen::SparseMatrix< real_t > positions_out;
    Eigen::SparseMatrix< real_t > du_out;
    Eigen::SparseMatrix< real_t > dv_out;
    compute_subdivision_coefficients_for_mesh(
        // Cube
        // CUBE,
        // Torus
         TORUS,
        us, vs,
        position_out, &du_out, &dv_out
        );

// 	std::cout << "position_out: " << position_out << "\n";
// 	std::cout << "du coefficients: " << du_out << "\n";
// 	std::cout << "dv coefficients: " << dv_out << "\n";

}

/*
void usage( const char* argv0 )
{
    std::cerr << "Usage: " << argv0 << " test_subdivision_matrices_eigen\n";
    exit( -1 );
}
int main( int argc, char* argv[] )
{
    using std::string;
    
    if( argc != 2 ) usage( argv[0] );
    
    const int resolution = atoi( argv[1] );
    if( resolution < 0 ) usage( argv[0] );
    
	Eigen::SparseMatrix< subdivision_matrix::real_t > position_coeffs;
	
    test_compute_coefficients( resolution, position_coeffs );

    return 0;
}

*/
