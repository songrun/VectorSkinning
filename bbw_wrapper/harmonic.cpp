/*
A "C" interface to libigl's harmonic coordinates.
*/

#include <cassert>

#define IGL_HEADER_ONLY
#include <igl/boundary_conditions.h>
#include <igl/normalize_row_sums.h>
#include <igl/harmonic.h>

#include <Eigen/Dense>

#include <iostream>

#define kVertexDimension 3

extern "C"
{

typedef double real_t;
typedef int index_t;

// Returns 0 for success, anything else is an error.
int harmonic(
    /// Input Parameters
    // 'vertices' is a pointer to num_vertices*kVertexDimension floating point values,
    // packed: x0, y0, z0, x1, y1, z1, ...
    // In other words, a num_vertices-by-kVertexDimension matrix packed row-major.
    int num_vertices, real_t* vertices,
    // 'faces' is a pointer to num_faces*3 integers,
    // where each face is three vertex indices: f0.v0, f0.v1, f0.v2, f1.v0, f1.v1, f1.v2, ...
    // Face i's vertices are: vertices[ faces[3*i]*2 ], vertices[ faces[3*i+1]*2 ], vertices[ faces[3*i+2]*2 ]
    // In other words, a num_faces-by-3 matrix packed row-major.
    int num_faces, index_t* faces,
    // 'boundary_indices' is a pointer to num_boundary_vertices integers,
    // where each element "i" in boundary_indices references the vertex whose data
    // is located at vertices[ boundary_indices[i]*kVertexDimension ].
    int num_boundary_indices, index_t* boundary_indices,
    // Power of the harmonic operation (1 is harmonic, 2 is bi-harmonic, etc ),
    int power,
    
    /// Output Parameters
    // 'Wout' is a pointer to num_vertices*num_boundary_indices values.
    // Upon return, W will be filled with each vertex in 'num_vertices' weight for
    // each boundary vertex in 'boundary_indices'.
    // The data layout is that all 'num_boundary_indices' weights for vertex 0
    // appear before all 'num_boundary_indices' weights for vertex 1, and so on.
    // In other words, a num_vertices-by-num_boundary_indices matrix packed row-major.
    real_t* Wout
    )
{
    using namespace std;
    using namespace igl;
    using namespace Eigen;
    
    assert( num_vertices > 0 );
    assert( vertices );
    assert( num_faces > 0 );
    assert( faces );
    assert( num_boundary_indices > 0 );
    assert( boundary_indices );
    
    assert( Wout );
    
    // #V by 2 list of mesh vertex positions
    // Make this an Eigen::map from 'vertices'
    MatrixXd V = Eigen::Map< Eigen::Matrix< real_t, Eigen::Dynamic, kVertexDimension, Eigen::RowMajor > >( vertices, num_vertices, kVertexDimension );
    // #F by 3 list of triangle indices
    // Make this an Eigen::map from 'faces'
    MatrixXi F = Eigen::Map< Eigen::Matrix< index_t, Eigen::Dynamic, 3, Eigen::RowMajor > >( faces, num_faces, 3 );
    
    // Make this an Eigen::map from boundary_indices.
    // VectorXi b = Eigen::Map< Eigen::Matrix< index_t, Eigen::Dynamic, 1 > >( boundary_indices, num_boundary_indices, 1 );
    // VectorXd bc( num_boundary_indices, num_boundary_indices );
    
    // Alec: If your mesh is (V,F) and your control cage is C then first build a #C by 2 MatrixXi of “edge indices” for your control cage. Easy, just:
    MatrixXd C(num_boundary_indices,kVertexDimension);
    MatrixXi CE(num_boundary_indices,2);
    VectorXi P(num_boundary_indices);
    for( int c = 0; c < num_boundary_indices; ++c )
    {
      C(c) = V( boundary_indices[c] );
      P(c) = c;
      CE(c,0) = c;
      CE(c,1) = (c+1)%C.rows();
    }
    // We don't use this, but we need an empty one to pass to boundary_conditions().
    MatrixXi BE;
    
    // Compute boundary conditions (aka fixed value constraints)
    // List of boundary indices (aka fixed value indices into VV)
    VectorXi b;
    // List of boundary conditions of each weight function
    MatrixXd bc;
    if(!boundary_conditions(V,F,C,P,BE,CE,b,bc))
    {
        return 1;
    }
    
    // compute harmonic coordinates
    // Weights matrix
    MatrixXd W;
    if(!harmonic(V,F,b,bc,power,W))
    {
        return 2;
    }
    
    // Normalize weights.
    normalize_row_sums(W,W);
    
    // Save output
    // Copy W to Wout.
    Eigen::Map< Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > >( Wout, num_vertices, num_boundary_indices ) = W;
    
    return 0;
}

}
