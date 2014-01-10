#include <cassert>
#include <vector>
#include <algorithm>

extern "C"
{
#include "raytri_wrapper.h"

int intersect_triangle1(double orig[3], double dir[3],
			double vert0[3], double vert1[3], double vert2[3],
			double *t, double *u, double *v);

struct raymesh_intersection_t* ray_mesh_intersections(
    double* ray_origin, double* ray_direction,
    int n_vertices, double* vertices, int n_faces, int* faces,
    int* num_intersections_out
    )
{
    assert( ray_origin );
    assert( ray_direction );
    assert( n_vertices >= 0 );
    assert( vertices );
    assert( n_faces >= 0 );
    assert( faces );
    assert( num_intersections_out );
    
    std::vector< raymesh_intersection_t > intersections;
    
    raymesh_intersection_t isect;
    int intersected = 0;
    for( int i = 0; i < n_faces; ++i )
    {
        intersected = intersect_triangle1(
            ray_origin, ray_direction,
            vertices + 3 * faces[ 3*i + 0 ],
            vertices + 3 * faces[ 3*i + 1 ],
            vertices + 3 * faces[ 3*i + 2 ],
            &isect.t,
            &isect.u,
            &isect.v
            );
        
        if( intersected )
        {
            isect.fi = i;
            intersections.push_back( isect );
        }
    }
    
    struct raymesh_intersection_t* result = new raymesh_intersection_t[ intersections.size() ];
    std::copy( intersections.begin(), intersections.end(), result );
    *num_intersections_out = intersections.size();
    return result;
}

void ray_mesh_intersections_free( struct raymesh_intersection_t* intersections )
{
    delete [] intersections;
}

}
