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

/* Our scalar and point type */
typedef double catmull_real_t;
typedef struct { catmull_real_t x, y, z; } coord_t, *coord;


/* A list type, used by model */
typedef struct {
	int n, alloc;
	void **buf;
} list_t, *list;

list list_new();
int list_push(list lst, void *e);
void list_del(list l);

#define len(l) ((l)->n)
#define elem(l, i) ((l)->buf[i])
#define foreach(i, e, l) for (i = 0, e = len(l) ? elem(l,i) : 0;	\
				i < len(l);				\
				e = (++i == len(l) ? 0 : elem(l, i)))


/* Geometry types */
typedef struct vertex_t {
	coord_t pos, avg_norm;
	list e, f;
	struct vertex_t * v_new;
	int idx;
} *vertex, vertex_t;
 
typedef struct {
	list f;
	vertex v[2];
	vertex e_pt; /* edge point for catmul */
	coord_t avg;
} *edge, edge_t;
 
typedef struct {
	list e, v;
	coord_t norm;
	vertex avg; /* average of all vertices, i.e. face point */
} *face, face_t;

/* The model type */
typedef struct { list e, v, f; } model_t, *model;

void model_norms(model m);
int model_add_vertex(model m, catmull_real_t x, catmull_real_t y, catmull_real_t z);
int model_add_face(model m, int n, ...);
int model_add_face_p(model m, int n, int* vs);
model model_new();
void model_del(model m);

/*
Performs one step of Catmull-Clark subdivision.
Returns a copy of 'm'.
*/
model catmull(model m);

model star();
model donut();
model cube();


#ifdef __cplusplus
}
#endif

#endif /* _CATMULL_H_ */
