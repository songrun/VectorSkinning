//------------------------------------------------------------------------------
//  Copyright 2014-2015 by Yotam Gingold and Songrun Liu
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

// #include "igl/readOBJ.h"
#include <iostream>
#include "catmull.h"
#include <stdio.h>
#include <vector>
#include <string.h>

void usage( char* argv0, int exit_code = -1 )
{
	std::cerr << "Usage: " << argv0 << "path/to/file.obj\n";
	exit( exit_code );
}

void error( std::string message, int exit_code = -1 )
{
	std::cerr << "Error: " << message << std::endl;
	exit( exit_code );
}

int cc_model_add_face(model m, std::vector<int > vs)
{
	int i, x;
	list lst = list_new();
 
	for (i = 0; i < vs.size(); i++) {

		list_push(lst, elem(m->v, vs[i]));
	}
 
	x = model_add_face_v(m, lst);
	list_del(lst);
	return x;
}

model readOBJ( const std::string obj_file_name )
{
	model m = model_new();
	
	FILE * obj_file = fopen(obj_file_name.c_str(),"r");
	if ( NULL==obj_file )
	{
		fprintf(stderr,"IOError: %s could not be opened...\n", obj_file_name.c_str());
		return m;
	}
	
	double x, y, z;
	char line[LINE_MAX];
	int line_no = 1;
	while (fgets(line, LINE_MAX, obj_file) != NULL) {
		char type[LINE_MAX];
		// Read first word containing type
	    if(sscanf(line, "%s",type) == 1) {
			// Get pointer to rest of line right after type
			char * l = &line[strlen(type)];
		
			if( strcmp(type, "v") == 0 )
			{
				int count = sscanf(l,"%lf %lf %lf\n",&x,&y,&z);
				if(count != 3)
				{
					fprintf(stderr, "Error: readOBJ() vertex on line %d should have 3 coordinates", 
					  line_no);
					fclose(obj_file);
					return m;
				}

				model_add_vertex(m, x, y, z);
			}
			else if( strcmp(type, "f") == 0 )
			{
				std::vector<int > f;
				std::vector<int > ftc;
				std::vector<int > fn;
				// Read each "word" after type
				char word[LINE_MAX];
				int offset;
				while(sscanf(l,"%s%n",word,&offset) == 1)
				{
					// adjust offset
					l += offset;
					// Process word
					unsigned int i,it,in;
					if (sscanf(word,"%u/%u/%u",&i,&it,&in) == 3)
					{
						f.push_back(i-1);
						ftc.push_back(it-1);
						fn.push_back(in-1);
					} else if(sscanf(word,"%u/%u",&i,&it) == 2)
					{
						f.push_back(i-1);
						ftc.push_back(it-1);
					} else if(sscanf(word,"%u//%u",&i,&in) == 2)
					{
						f.push_back(i-1);
						fn.push_back(in-1);
					} else if(sscanf(word,"%u",&i) == 1)
					{
						f.push_back(i-1);
					} else
					{
						fprintf(stderr,
							"Error: readOBJ() face on line %d has invalid element format\n",
							line_no);
						fclose(obj_file);
						return m;
					}
				}
				if((f.size()>0 && fn.size() == 0 && ftc.size() == 0) ||
				   (f.size()>0 && fn.size() == f.size() && ftc.size() == 0) ||
				   (f.size()>0 && fn.size() == 0 && ftc.size() == f.size()) ||
				   (f.size()>0 && fn.size() == f.size() && ftc.size() == f.size()))
				{
					// No matter what add each type to lists so that lists are the
					// correct lengths
					cc_model_add_face(m, f);
// 					FTC.push_back(ftc);
// 					FN.push_back(fn);
				} 
				else {
					fprintf(stderr,
						  "Error: readOBJ() face on line %d has invalid format\n", line_no);
					fclose(obj_file);
					return m;
				}
			} else if(strlen(type) >= 1 && (type[0] == '#' || 
				type[0] == 'g'  ||
				type[0] == 's'  ||
				strcmp("usemtl",type)==0 ||
				strcmp("mtllib",type)==0))
			{
			//ignore comments or other shit
			}
		}
		line_no++;
	}
	fclose(obj_file);
	
	model_norms(m);
	return m;
	
}


int main(int argc, char **argv)
{
	int i;
	void *p;

	models = list_new();
	
	model obj = readOBJ( argv[1] );
	list_push(models, obj);
	model_pos = 1;

	init_gfx(&argc, argv);

	foreach(i, p, models) {
		model model_p = reinterpret_cast<model>(p);
		model_del(model_p);
	}
	list_del(models);

	return 0;
	
}
/*
int main( int argc, char ** argv ) {
	
	if( 2 != argc )
	{
		usage(argv[0]);
	}
	
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	const int success = igl::readOBJ( argv[1], V, F );
	if( !success ) usage(argv[0]);
	
// 	cout << "V:\n" << V << '\n';
// 	cout << "F:\n" << F << '\n';
	
	Model G;
	bool is_success = G.build( V, F );
	if( !is_success ) error("Fail to build model.");
	
 	G.print();
 	
// 	is_success = G.subdivide(2);
// 	if( !is_success ) error("Fail to subdivide model.");
// 	
// 	G.print();
// 	
	return 0;
}
*/