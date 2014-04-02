#include "catmull_eigen.h"

#include <cassert>

void model_faces_to_eigen(
	model m,
	Eigen::MatrixNi& faces_out,
	Eigen::MatrixNd& face_normals_out
	)
{
	faces_out.set_size( len(m->f), 4 );
	face_normals_out.set_size( len(m->f), 3 );
	
	int face_index, face_vertex_index;
	void *rf, *rv;
	face f;
	vertex v;
	foreach(face_index, rf, m->f) {
		f = static_cast< face >( rf );
		
		face_normals_out[ face_index, 0 ] = f->norm.x;
		face_normals_out[ face_index, 1 ] = f->norm.y;
		face_normals_out[ face_index, 2 ] = f->norm.z;
		
		// This code assumes quad faces.
		assert( 4 == len(f->v) );
		
		foreach(face_vertex_index, rv, f->v) {
			v = static_cast< vertex >( rv );
			
			faces_out[ face_index, face_vertex_index ] = v->idx;
		}
	}
}

void model_vertices_to_eigen(
	model m,
	Eigen::MatrixNd& vertex_positions_out,
	Eigen::MatrixNd& vertex_normals_out
	)
{
	vertex_positions_out.set_size( len(m->v), 3 );
	vertex_normals_out.set_size( len(m->v), 3 );
	
	int vertex_index;
	void *rv;
	vertex v;
	
	foreach(vertex_index, rv, m->v) {
		v = static_cast< vertex >( rv );
		
		vertex_positions_out[ vertex_index, 0 ] = v->pos.x;
		vertex_positions_out[ vertex_index, 1 ] = v->pos.y;
		vertex_positions_out[ vertex_index, 2 ] = v->pos.z;
		
		vertex_normals_out[ vertex_index, 0 ] = v->avg_norm.x;
		vertex_normals_out[ vertex_index, 1 ] = v->avg_norm.y;
		vertex_normals_out[ vertex_index, 2 ] = v->avg_norm.z;
	}
}
			vertex_normals_out[ 
			if (interp_norm)
				glNormal3dv(&(v->avg_norm.x));
			glVertex3dv(&(v->pos.x));
		}
		glEnd();
	}
}
