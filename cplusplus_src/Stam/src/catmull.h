//------------------------------------------------------------------------------
//  From: http://rosettacode.org/wiki/Catmull%E2%80%93Clark_subdivision_surface/C
//  Modified 2014-2015 by Yotam Gingold and Songrun Liu
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#ifndef _CATMULL_H_
#define _CATMULL_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

int max_depth = 7;
int show_parent = 1;
int wireframe_mode = 1;
int face_mode = 1;
int model_idx = 0;
int interp_norm = 0;
int LINE_MAX = 100000;

#define new_struct(type, var) type var = calloc(sizeof(type##_t), 1)
#define len(l) ((l)->n)
#define elem(l, i) ((l)->buf[i])
#define foreach(i, e, l) for (i = 0, e = len(l) ? elem(l,i) : 0;	\
				i < len(l);				\
				e = (++i == len(l) ? 0 : elem(l, i)))

typedef struct {
	int n, alloc;
	void **buf;
} list_t, *list;

typedef struct { list e, v, f; } model_t, *model;

list list_new();
int list_push(list lst, void *e);
void init_gfx(int *c, char **v) ;
void list_del(list l);
void model_del(model m);

model cube();

// global
list models = 0;
int model_pos = 0;

void model_norms(model m);
int model_add_vertex(model m, GLfloat x, GLfloat y, GLfloat z);
int model_add_face_v(model m, list vl);
model model_new();

#ifdef __cplusplus
}
#endif

#endif /* _CATMULL_H_ */