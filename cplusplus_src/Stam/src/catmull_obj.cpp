#include "catmull_obj.h"

#include <cstdio>
#include <vector>
#include <iostream>
#include <cstring>

#define LINE_MAX 8192

model readOBJ( const char* obj_file_name )
{
	model m = model_new();
	
	FILE * obj_file = fopen(obj_file_name,"r");
	if ( NULL==obj_file )
	{
		fprintf(stderr,"IOError: %s could not be opened...\n", obj_file_name);
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
					model_add_face_p(m, f.size(), &f[0]);
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
