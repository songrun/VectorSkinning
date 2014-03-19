//------------------------------------------------------------------------------
//  Copyright 2014-2015 by Yotam Gingold and Songrun Liu
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#include "Model.h"

bool Model::build( const Eigen::MatrixXd V, const Eigen::MatrixXi F )
{
	uint N = F.cols();
	
	// allocate memory
	v_size = V.rows();
	f_size = F.rows();
	e_size = F.size()/2; 	
	vertices = new Vertex[v_size];   
    faces = new Face[f_size];
    edges = new Edge[e_size];
    assert ( vertices && faces && edges );
    
    for ( uint i=0; i<v_size; i++ ) {
    	vertices[i].p = V.row(i);
    }
	
	uint e_count = 0;
	for ( uint i=0; i<f_size; i++ ) {
		
		uint v_id = F( i, N-1 );
		for ( uint j=0; j<N; j++ ) {			
			// search for existence of current edge
			uint e_id = 0;
			for ( uint k=0; k<e_count; k++ ) {
				// the edge does exist
				if ( (edges[k].vid[0] == v_id && edges[k].vid[1] == F(i,j)) ||
					 (edges[k].vid[1] == v_id && edges[k].vid[0] == F(i,j)) ) {
					e_id = k;
					break; 
				}
				// otherwise
				else
					e_id++;
			}
			
			edges[e_id].vid[0] = v_id;
			edges[e_id].vid[1] = F(i,j);
			
			// one of the fid should be UINT_MAX before the assignment
			assert ( UINT_MAX == edges[e_id].fid[0] || UINT_MAX == edges[e_id].fid[1] );
			if ( UINT_MAX == edges[e_id].fid[0] )	
				edges[e_id].fid[0] = i;
			else
				edges[e_id].fid[1] = i;

			// edge vector
// 			edges[e_id].v = vertices[edges[e_id].vid[1]].p - vertices[edges[e_id].vid[0]].p; 
			// end of edge assignment
			
			// if current edge is a new edge
			if ( e_id == e_count )	{
				e_count++;
				vertices[v_id].v_e.push_back( e_id );
				vertices[F(i,j)].v_e.push_back( e_id );
			}	
			
			faces[i].es.push_back( e_id );

			v_id = F(i,j);
			faces[i].vs.push_back( v_id );
			vertices[v_id].v_f.push_back( i );
			
		}
	}

	return true;
}

void Model::print()
{
	cout << "print vertices: " << endl;
	for ( uint vi=0; vi<v_size; vi++ ) {
		cout << vertices[vi].p << endl;
	}
	
	cout << "print faces: " << endl;
	for ( uint fi=0; fi<f_size; fi++ ) {
		Face f = faces[fi];
// 		for ( list<uint>::iterator f_v=f.vs.begin(); f_v!=f.vs.end(); f_v++ ) {
//			cout << *f_v << ' ';
		for ( auto f_v : f.vs ) {
			cout << f_v << ' ';
		} 
		cout << endl;
	}
}


bool Model::subdive( const uint times )
{
	for ( uint i=0; uint<times; uint++ ) {
		// allocate space for subdivided model
		Vertex * new_vertices = new Vertex[v_size+f_size+e_size];
		uint count = 0;
		for	( uint i=0; i<f_size; i++ ) {
			assert( faces[i].vs.size() == faces[i].es.size() ) ;
			count += faces[i].vs.size();
		}
		Face * new_faces = new Face[count];
		Edge * new_edges = new Edge[count+e_size*2];
		assert( new_vertices && new_faces && new_edges); 
		
		// Set each face point to be the average of all original points for the respective face.
		for	( uint i=0; i<f_size; i++ ) {
			Face f = faces[i];
			for ( list<uint>::iterator f_v=f.vs.begin(); f_v!=f.vs.end(); f_v++ ) {
				uint v_id = *f_v;
				new_vertices[v_size+i].p += vertices[v_id].p;
			}
			new_vertices[v_size+i].p /= f.vs.size();
		}
		
		// Set each edge point to be the average of the two neighbouring face points and its two original endpoints.
		for	( uint i=0; i<e_size; i++ ) {
			Edge e = edges[i];
			new_vertices[i+v_size+e_size].p = (e.vid[0].p + e.vid[1].p 
								+ f_vertices[e.fid[0]].p + f_vertices[e.fid[1]].p)/4;
			// splited edges
			new_edges[i*2].vid[0] = edges[i].vid[0];
			new_edges[i*2+1].vid[1] = edges[i].vid[1];
			new_edges[i*2].vid[1] = new_edges[i*2+1].vid[0] = v_size+f_size+i;
		}
		
		// For each face point, add an edge for every edge of the face, 
		// connecting the face point to each edge point for the face.
		uint f_count = 0;
		for ( uint i=0; i<f_size; i++ ) {
			Face f = faces[i];
			uint num = f.vs.size();
			uint e_count = f_count;
			
			for ( auto f_v : f.vs ) {
				new_faces[f_count].vs.push_back(v_size+i);
				new_faces[f_count].vs.push_back(f_v);
				
				for ( auto f_e : f.es ) {
					if ( f_v == edges[f_e].vid[0]  ) {
						new_faces[f_count].es.push_back(f_e*2);
						new_faces[f_count].vs.push_back(new_edges[2*f_e].vid[1]);
						new_vertices[new_edges[2*f_e].vid[1]] 
						
						// connect the face point to one edge point
						new_edges[2*count+e_count+j].vid[0] = new_edges[2*f_e].vid[1];
						new_edges[2*count+e_count+j].vid[1] = v_size+i;
						new_edges[2*count+e_count+j].fid[0] = f_count;
						
					}
					else if ( f_v == edges[f_e].vid[1] ) {
						new_faces[f_count].es.push_back(f_e*2+1);
						new_faces[f_count].vs.push_back(new_edges[2*f_e+1].vid[0]);
						
					}
				}
				f_count++;
			}
				
		}
		
	}
	
	return true;
}
*/