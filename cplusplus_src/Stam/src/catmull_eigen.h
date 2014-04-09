#ifndef _CATMULL_EIGEN_H_
#define _CATMULL_EIGEN_H_

#include <Eigen/Core>

// Erases faces_out and fills it with the 4 indices of the faces of the model.
// Erases faces_normals_out and fills it with the xyz coordinates of each face's normal.
// NOTE: The model must have only quad faces for this to work.
void model_quad_faces_to_eigen(
	model m,
	Eigen::MatrixNi& faces_out,
	Eigen::MatrixNd& face_normals_out
	);

// Erases vertex_positions_out and vertex_normals_out and fills them
// with the 3 xyz coordinates of the vertices of the model and their normals.
void model_vertices_to_eigen(
	model m,
	Eigen::MatrixNd& vertex_positions_out,
	Eigen::MatrixNd& vertex_normals_out
	);

// Returns the indices of the control points of the given face.
std::vector< int > control_point_indices_of_face( model m, int face_index );

#endif /* _CATMULL_EIGEN_H_ */
