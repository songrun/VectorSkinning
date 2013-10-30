/*
A "C" interface to libigl's bbw.
*/

#include <cassert>

#define IGL_HEADER_ONLY
#define IGL_NO_MOSEK
#include <igl/boundary_conditions.h>
#include <igl/bbw/bbw.h>

#include <Eigen/Dense>

#include <iostream>

extern "C"
{

typedef double real_t;

// Returns 0 for success, anything else is an error.
int bbw2D(
    /// Input Parameters
    // 'vertices' is a pointer to num_vertices*2 floating point values,
    // packed: x0, y0, x1, y1, ...
    int num_vertices, real_t* vertices,
    // 'faces' is a pointer to num_faces*3 integers,
    // where each face is three vertex indices: f0.v0, f0.v1, f0.v2, f1.v0, f1.v1, f1.v2, ...
    // Face i's vertices are: vertices[ faces[3*i]*2 ], vertices[ faces[3*i+1]*2 ], vertices[ faces[3*i+2]*2 ]
    int num_faces, int* faces,
    // 'skeleton_vertices' is a pointer to num_skeleton_vertices*2 floating point values,
    // packed the same way as 'vertices'
    int num_skeleton_vertices, real_t* skeleton_vertices,
    // 'skeleton_point_handles' is a pointer to num_skeleton_point_handles integers,
    // where each element "i" in skeleton_point_handles references the vertex whose data
    // is located at skeleton_vertices[ skeleton_point_handles[i]*2 ].
    int num_skeleton_point_handles, int* skeleton_point_handles,
    // TODO: Take skeleton bone edges and cage edges
    
    /// Output Parameters
    // 'Wout' is a pointer to num_vertices*num_skeleton_vertices values.
    // Upon return, W will be filled with each vertex in 'num_vertices' weight for
    // each skeleton vertex in 'num_skeleton_vertices'.
    // The data layout is that all num_vertices values for skeleton vertex 0
    // appear before all num_vertices values for skeleton vertex 1, and so on.
    // TODO: This data layout may be transposed, please check.
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
    assert( num_skeleton_vertices > 0 );
    assert( skeleton_vertices );
    
    // TODO: These next two asserts need not be true once we add support for the other
    //       kinds of handles.
    assert( num_skeleton_point_handles > 0 );
    assert( skeleton_point_handles );
    
    assert( Wout );
    
    // #V by 3 list of mesh vertex positions
    // TODO: Make this an Eigen::map from 'vertices'
    MatrixXd V;
    // #F by 3 list of triangle indices
    // TODO: Make this an Eigen::map from 'faces'
    MatrixXi F;
    
    // "Skeleton" (handles) descriptors:
    // List of control and joint (bone endpoint) positions
    // TODO: Make this an Eigen::map from skeleton_vertices
    MatrixXd C;
    // List of point handles indexing C
    // TODO: Make this an Eigen::map from skeleton_point_handles
    VectorXi P;
    // List of bone edges indexing C
    // TODO: ...
    MatrixXi BE;
    // List of cage edges indexing *P*
    // TODO: ...
    MatrixXi CE;
    
    // Compute boundary conditions (aka fixed value constraints)
    // List of boundary indices (aka fixed value indices into VV)
    VectorXi b;
    // List of boundary conditions of each weight function
    MatrixXd bc;
    if(!boundary_conditions(V,F,C,P,BE,CE,b,bc))
    {
        return 1;
    }
    
    cout<<"b=["<<b<<"];"<<endl;
    cout<<"bc=["<<bc<<"];"<<endl;
    
    // compute BBW 
    // Default bbw data and flags
    BBWData bbw_data;
    // Weights matrix
    MatrixXd W;
    if(!bbw(V,F,b,bc,bbw_data,W))
    {
        return 1;
    }
    
    // Save output
    // TODO: Copy W to Wout.
    
    return 0;
}

}
