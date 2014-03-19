//------------------------------------------------------------------------------
//  Copyright 2014-2015 by Yotam Gingold and Songrun Liu
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#ifndef _SR_MODEL_H_
#define _SR_MODEL_H_

#include <list>
#include <Eigen/Dense>
#include <cassert>
#include <string>
#include <iostream>

typedef unsigned int uint;
typedef Eigen::RowVector3d Point3d;
typedef Eigen::RowVector3d Vector3d;

using namespace std;

struct Face
{
	Face(){ }
	
	list<uint> vs; // vertex id
 	list<uint> es; // edge id 
 	
};

struct Vertex
{
	Vertex() { p << 0,0,0; }
	
	Point3d p; // position

	list<uint> v_e; //a list of edges
    list<uint> v_f; //a list of faces
};

struct Edge
{
	Edge() { 
	vid[0]=vid[1]=UINT_MAX; 
	fid[0]=fid[1]=UINT_MAX;
	}
	
	uint vid[2];
	uint fid[2];
	
// 	Vector3d v; // parallel vector
};

struct Model
{
	Model() {
		v_size = e_size = f_size = 0;
        vertices=NULL;
        edges=NULL;
        faces=NULL; 
	}
	
	~Model() {
		delete [] faces;    faces=NULL;
        delete [] edges;    edges=NULL;
        delete [] vertices; vertices=NULL;
        v_size=e_size=f_size=0;
	}
	
	//build this model
    bool build( const Eigen::MatrixXd V, const Eigen::MatrixXi F );
    bool subdive( const uint times );
    void print();
    
    Vertex * vertices;  //vertices
    Face * faces;      //triangles
    Edge * edges;     //edges
    uint v_size, e_size, f_size;
};


#endif	// _SR_MODEL_H_