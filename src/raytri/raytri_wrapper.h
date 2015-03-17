#ifndef __raytri_wrapper_h__
#define __raytri_wrapper_h__

// Returns a pointer to a list of intersections between the ray emanating from the 3d point
// 'ray_origin' in the direction of the 3d vector 'ray_direction' with the mesh defined
// by the flattened list of vertices (succession of 'n_vertices' xyz coordinates)
// and the flattened list of faces (succession of 'n_faces' vertex indices).
// The number of intersections is returned in the integer output parameter 'num_intersections_out'.
// Each intersection has a 't' along the ray, the index 'fi' of the face in 'faces',
// and coordinates inside the triangle 'u' and 'v'.
struct raymesh_intersection_t
{
    double t;
    int fi;
    double u;
    double v;
};
struct raymesh_intersection_t* ray_mesh_intersections(
    double* ray_origin, double* ray_direction,
    int n_vertices, double* vertices, int n_faces, int* faces,
    int* num_intersections_out
    );
// Call this to free the memory returned by ray_mesh_intersections()
void ray_mesh_intersections_free( struct raymesh_intersection_t* intersections );

#endif /* __raytri_wrapper_h__ */
