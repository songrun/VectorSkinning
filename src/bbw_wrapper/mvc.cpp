/*
A "C" interface to libigl's mvc.
*/

#include <cassert>

#define IGL_HEADER_ONLY
#include <igl/mvc.h>

#include <iostream>

#define kVertexDimension 2

extern "C"
{

typedef double real_t;

// Returns 0 for success, anything else is an error.
int mvc(
    /// Input Parameters
    // 'vertices' is a pointer to num_vertices*2 floating point values,
    // packed: x0, y0, x1, y1, ...
    // In other words, a num_vertices-by-2 matrix packed row-major.
    int num_vertices, real_t* vertices,
    // 'line_loop' is a pointer to num_line_loop*2 floating point values,
    // packed: x0, y0, x1, y1, ...
    // In other words, a num_line_loop-by-2 matrix packed row-major.
    int num_line_loop, real_t* line_loop,
    
    /// Output Parameters
    // 'Wout' is a pointer to num_vertices*num_line_loop values.
    // Upon return, W will be filled with each vertex in 'num_vertices' weight for
    // each vertex in 'line_loop'.
    // The data layout is that all 'num_line_loop' weights for vertex 0
    // appear before all 'num_line_loop' weights for vertex 1, and so on.
    // In other words, a num_vertices-by-num_line_loop matrix packed row-major.
    real_t* Wout
    )
{
    using namespace std;
    using namespace Eigen;
    
    assert( num_vertices > 0 );
    assert( vertices );
    assert( num_line_loop >= 3 );
    assert( line_loop );
    
    assert( Wout );
    
    // #V by 2 list of mesh vertex positions
    // Make this an Eigen::map from 'vertices'
    MatrixXd V = Eigen::Map< Eigen::Matrix< real_t, Eigen::Dynamic, kVertexDimension, Eigen::RowMajor > >( vertices, num_vertices, kVertexDimension );
    // #H by 2 list of cage vertex positions
    // Make this an Eigen::map from 'line_loop'
    MatrixXd C = Eigen::Map< Eigen::Matrix< real_t, Eigen::Dynamic, kVertexDimension, Eigen::RowMajor > >( line_loop, num_line_loop, kVertexDimension );
    
    MatrixXd W;
    igl::mvc(V,C,W);
    
    // Normalize weights.
    // Q: Do we need this?
    // A: No. igl::mvc() does it.
    // normalize_row_sums(W,W);
    
    // Save output
    // Copy W to Wout.
    Eigen::Map< Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > >( Wout, num_vertices, num_line_loop ) = W;
    
    return 0;
}

}
