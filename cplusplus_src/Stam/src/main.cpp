//------------------------------------------------------------------------------
//  Copyright 2014-2015 by Yotam Gingold and Songrun Liu
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

// #include "igl/readOBJ.h"
#include <iostream>
#include "catmull.h"
#include "catmull_obj.h"
#include "stam_eigen_reader.h"
#include <stdio.h>
#include <vector>
#include <string.h>

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

namespace
{
int gwin;
GLfloat rot[] = { 20, 40, 0 };
GLfloat rotx = 0, roty = 1;
GLfloat ambient[] = { .5f, .5f, .5f, 1.f };
GLfloat diffuse[] = { .5f, .5f, .5f, .6f };
GLfloat litepos[] = { 0, 2, 3, 1 };
GLfloat litepos2[] = { 0, -2, 5, 1 };
GLfloat color[] = { .3, 1, .6 };
GLfloat hole[] = { .7, .4, 0 };
GLfloat color2[] = { .5, .5, .5 };
GLfloat hole2[] = { .5, .2, 0 };
GLfloat red[] = {1, 0, 0};
GLfloat zpos = -6;
 
int max_depth = 7;
int show_parent = 1;
int wireframe_mode = 1;
// int face_mode = 1;
int model_idx = 0;
int interp_norm = 0;

list models = 0;
int model_pos = 0;
}

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

void draw_model(model m)
{
	int i, j;
	void *rf, *rv;
	face f;
	vertex v;
	foreach(i, rf, m->f) {
		f = static_cast< face >( rf );
		glBegin(GL_POLYGON);
		if (!interp_norm) glNormal3dv(&(f->norm.x));
		foreach(j, rv, f->v) {
			v = static_cast< vertex >( rv );
			if (interp_norm)
				glNormal3dv(&(v->avg_norm.x));
			glVertex3dv(&(v->pos.x));
		}
		glEnd();
	}
}
 
void draw_wireframe(model m, GLfloat *color, GLfloat *hole_color)
{
	int i;
	edge e;
	void *re;
 
	glDisable(GL_LIGHTING);
	foreach(i, re, m->e) {
		e = static_cast< edge >( re );
		if (e->f->n != 2) glColor3fv(hole_color);
		else		  glColor3fv(color);
 
		glBegin(GL_LINES);
		glVertex3dv(&(e->v[0]->pos.x));
		glVertex3dv(&(e->v[1]->pos.x));
		glEnd();
	}
}
 
void draw_faces(model m)
{
	glPushMatrix();
	glLoadIdentity();
	glEnable(GL_LIGHTING);
	glLightfv(GL_LIGHT0, GL_AMBIENT,  ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE,  diffuse);
	glLightfv(GL_LIGHT0, GL_POSITION, litepos);
	glEnable(GL_LIGHT0);
 
	glLightfv(GL_LIGHT1, GL_DIFFUSE,  diffuse);
	glLightfv(GL_LIGHT1, GL_POSITION, litepos2);
	glEnable(GL_LIGHT1);
	glPopMatrix();
 
	if (wireframe_mode) {
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.0, 1.0);
	}
	draw_model(m);
	if (wireframe_mode)
		glDisable(GL_POLYGON_OFFSET_FILL);
}
 
void keyspecial(int k, int x, int y)
{
	switch(k) {
	case GLUT_KEY_UP:
		rotx --; break;
	case GLUT_KEY_DOWN:
		rotx ++; break;
	case GLUT_KEY_LEFT:
		roty --; break;
	case GLUT_KEY_RIGHT:
		roty ++; break;
	}
}
 
void set_model()
{
	int i;
	void *p;
 
	model_pos = 1;
	model_idx = (model_idx + 1) % 3;
 
	foreach(i, p, models) model_del( static_cast< model >( p ) );
 
	len(models) = 0;
 
	switch(model_idx) {
	case 0:
		list_push(models, cube()); break;
	case 1:
		list_push(models, donut()); break;
	case 2:
		list_push(models, star()); break;
	}
}
 
void keypress(unsigned char key, int x, int y)
{
	switch(key) {
	case 27: case 'q':
		glFinish();
		glutDestroyWindow(gwin);
		return;
	case 'w':
		wireframe_mode = (wireframe_mode + 1) % 3;
		return;
	case 'l':
		diffuse[0] += .1;
		diffuse[1] += .1;
		diffuse[2] += .1;
		return;
	case 'L':
		diffuse[0] -= .1;
		diffuse[1] -= .1;
		diffuse[2] -= .1;
		return;
	case ' ':
		rotx = roty = 0;
		return;
	case 'z':
		zpos ++;
		return;
	case 'a':
		zpos --;
		return;
	case 's':
		interp_norm = !interp_norm;
		break;
	case 'p':
		show_parent = (show_parent + 1) % 3;
		break;
 
	case '.': case '>':
		if (++model_pos >= max_depth) model_pos = max_depth;
		return;
	case ',': case '<':
		if (--model_pos < 0) model_pos = 0;
		return;
 
	case 'm':
		set_model();
		break;
	}
}
 
void render()
{
	if (!len(models)) return;
	while (model_pos >= len(models))
		list_push(models, catmull(static_cast< model >( elem(models, len(models) - 1))));
 
	glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
 
	glTranslatef(.0f, 0.f, zpos);
 
	rot[0] += rotx / 8.;
	rot[1] += roty / 8.;
 
	if (rot[0] > 360) rot[0] -= 360;
	if (rot[0] < 0) rot[0] += 360;
	if (rot[1] > 360) rot[1] -= 360;
	if (rot[1] < 0) rot[1] += 360;
 
	glRotatef(rot[0], 1, 0, 0);
	glRotatef(rot[1], 0, 1, 0);
 
	GLfloat *pcolor = color2;
	if (model_pos && show_parent) {
		if (show_parent == 2) pcolor = red;
		draw_wireframe(static_cast< model >(elem(models, model_pos - 1)), pcolor, hole2);
	}
 
	model m = static_cast< model >( elem(models, model_pos) );
	if (wireframe_mode) draw_faces(m);
	if (wireframe_mode < 2) draw_wireframe(m, color, hole);
 
	glFlush();
	glFinish();
	glutSwapBuffers();
}
 
void resize(int w, int h)
{
	printf("size %d %d\n", w, h);
	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45.f, (GLfloat)w / h, .1f, 100.f);
	glMatrixMode(GL_MODELVIEW);
}
 
void init_gfx(int *c, char **v) {
	glutInit(c, v);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);
	glutInitWindowSize(640, 480);
 
	gwin = glutCreateWindow("Catmull-Clark");
	glutReshapeFunc(resize);
 
	glClearColor(.0f, .0f, .0f, .0f);
	glClearDepth(1.0f);
	glDepthFunc(GL_LESS);
	glEnable(GL_DEPTH_TEST);
	glShadeModel(GL_SMOOTH);
 
	glutKeyboardFunc(keypress);
	glutSpecialFunc(keyspecial);
	glutDisplayFunc(render);
	glutIdleFunc(render);
	glutMainLoop();
}

int main(int argc, char **argv)
{
	EVALSTRUCT ** ev;
	int Nmax;

	ev = read_eval ( &Nmax );
	if ( !ev ) {
		printf ("Error while reading from Stam data.\n"); 
		exit ( 1 );
	}

	print_eval ( ev, Nmax );
	exit ( 0 );

	return 0;
}
   
 /*  
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
		model model_p = static_cast<model>(p);
		model_del(model_p);
	}
	list_del(models);

	return 0;
	
}

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
