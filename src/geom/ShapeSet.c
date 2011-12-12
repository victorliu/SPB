#include "ShapeSet.h"
#include "BVH.h"
#include <stdlib.h>
#include <string.h>

#define ShapeSet2_threshold 32
#define ShapeSet3_threshold 32

//#define SS_DEBUG

#ifdef SS_DEBUG
# include <stdio.h>
# define SSDBG(...) fprintf(stderr, __VA_ARGS__)
#else
# define SSDBG(...)
#endif


typedef struct shape2_list_{
	shape2 s;
	int tag;
	struct shape2_list_ *next;
} *shape2_list;


typedef struct{
	shape2 *s;
	aabb2 b;
	int tag;
} vecelem2;

struct ShapeSet2_{
	unsigned int n;
	shape2_list lst;
	
	vecelem2 *v; // always exists
	BVH2 bvh; // may not exist
	
	int finalized;
	int use_bvh;
	
	int periodic;
	double lattice[4];
};

ShapeSet2 ShapeSet2_new(const double lattice[4]){
	ShapeSet2 ss = (ShapeSet2)malloc(sizeof(struct ShapeSet2_));
	ss->n = 0;
	ss->lst = NULL;
	ss->bvh = NULL;
	ss->v = NULL;
	ss->finalized = 0;
	ss->use_bvh = 0;
	ss->periodic = 0;
	if(NULL != lattice){
		ss->periodic = 1;
		ss->lattice[0] = lattice[0];
		ss->lattice[1] = lattice[1];
		ss->lattice[2] = lattice[2];
		ss->lattice[3] = lattice[3];
	}
	return ss;
}

static void ShapeSet2_unfinalize(ShapeSet2 ss){
	if(ss->finalized){
		if(ss->use_bvh){
			BVH2_destroy(ss->bvh);
		}else{
			free(ss->v);
		}
		ss->finalized = 0;
	}
}

void ShapeSet2_destroy(ShapeSet2 ss){
	if(NULL == ss){ return; }
	ShapeSet2_unfinalize(ss);
	while(NULL != ss->lst){ // free the linked list
		shape2_list t = ss->lst;
		ss->lst = ss->lst->next;
		free(t);
	}
	free(ss);
}

int ShapeSet2_add(ShapeSet2 ss, const shape2 *s, int tag){
	unsigned int i, ilim = 1;
	static const int off[] = {
		0, 0, // must be first
		 1,  0,
		-1,  0,
		 0,  1,
		 0, -1,
		 1,  1,
		 1, -1,
		-1,  1,
		-1, -1
	};
	if(NULL == ss){ return -1; }
	if(NULL == s){ return -2; }
	
	ShapeSet2_unfinalize(ss);
	
	if(ss->periodic){
		ilim += 8;
	}
	for(i = 0; i < ilim; ++i){
		shape2_list t = (shape2_list)malloc(sizeof(struct shape2_list_));
		memcpy(&(t->s), s, sizeof(shape2));
		SSDBG("Adding shape type=%d\n", t->s.type);
		t->s.o[0] += off[2*i+0] * ss->lattice[0] + off[2*i+1] * ss->lattice[2];
		t->s.o[1] += off[2*i+0] * ss->lattice[1] + off[2*i+1] * ss->lattice[3];
		t->tag = tag;
		t->next = ss->lst;
		ss->lst = t;
		ss->n++;
	}
}


typedef struct{
	int index;
	const vecelem2 *v;
} shape2_iter_data;
static int shape2_iter(double b[4], int *tag, void *data){
	shape2_iter_data *d = (shape2_iter_data*)data;
	*tag = d->index; // the BVH's tag is actuall the index into the vector v
	b[0] = d->v[d->index].b.c[0] - d->v[d->index].b.h[0];
	b[1] = d->v[d->index].b.c[0] + d->v[d->index].b.h[0];
	b[2] = d->v[d->index].b.c[1] - d->v[d->index].b.h[1];
	b[3] = d->v[d->index].b.c[1] + d->v[d->index].b.h[1];
	d->index++;
	return 1;
}
static int vecelem2_cmp(const void *a, const void *b){
	const vecelem2 *e0 = (const vecelem2*)a;
	const vecelem2 *e1 = (const vecelem2*)b;
	return e1->tag - e0->tag;
}

static void ShapeSet2_finalize(ShapeSet2 ss){
	if(ss->finalized){ return; }
	
	ss->v = (vecelem2*)malloc(sizeof(vecelem2)*ss->n);
	// Insert linked list elements into vector v
	{
		shape2_list t = ss->lst;
		int index = 0;
		while(NULL != t){
			ss->v[index].s = &(t->s);
			shape2_get_aabb(&(t->s), &(ss->v[index].b));
			ss->v[index].tag = t->tag;
			++index;
			t = t->next;
		}
	}
	// Sort vector v so largest tags are at lower offsets
	qsort(ss->v, ss->n, sizeof(vecelem2), &vecelem2_cmp);
	
	if(ss->n >= ShapeSet2_threshold){
		shape2_iter_data d;
		d.index = 0;
		d.v = ss->v;
		ss->bvh = BVH2_new(ss->n, &shape2_iter, (void*)&d);
		ss->use_bvh = 1;
	}
	ss->finalized = 1;
}

static int query_pt2(int tag, const double b[4], void *data){
	int *best = (int*)data;
	if(tag > *best){ *best = tag; }
}

int ShapeSet2_query_pt(ShapeSet2 ss, const double p[2], const shape2 **s, int *tag){
	if(NULL == ss){ return -1; }
	if(NULL == p){ return -2; }
	
	ShapeSet2_finalize(ss);
	
	if(ss->use_bvh){
		int best = -1;
		BVH2_query_pt(ss->bvh, p, &query_pt2, &best);
		if(best >= 0){
			if(shape2_contains(ss->v[best].s, p)){
				if(NULL != tag){
					*tag = ss->v[best].tag;
				}
				if(NULL != s){
					*s = ss->v[best].s;
				}
				return 1;
			}
		}
	}else{
		int i;
		for(i = 0; i < ss->n; ++i){
			if(aabb2_contains(&(ss->v[i].b), p)){
				if(shape2_contains(ss->v[i].s, p)){
					if(NULL != tag){
						*tag = ss->v[i].tag;
					}
					if(NULL != s){
						*s = ss->v[i].s;
					}
					return 1;
				}
			}
		}
	}
	return 0;
}

























typedef struct shape3_list_{
	shape3 s;
	int tag;
	struct shape3_list_ *next;
} *shape3_list;


typedef struct{
	shape3 *s;
	aabb3 b;
	int tag;
} vecelem3;

struct ShapeSet3_{
	unsigned int n;
	shape3_list lst;
	
	vecelem3 *v; // always exists
	BVH3 bvh; // may not exist
	
	int finalized;
	int use_bvh;
	
	int periodic;
	double lattice[9];
};

ShapeSet3 ShapeSet3_new(const double lattice[9]){
	ShapeSet3 ss = (ShapeSet3)malloc(sizeof(struct ShapeSet3_));
	ss->n = 0;
	ss->lst = NULL;
	ss->bvh = NULL;
	ss->v = NULL;
	ss->finalized = 0;
	ss->use_bvh = 0;
	ss->periodic = 0;
	if(NULL != lattice){
		ss->periodic = 1;
		ss->lattice[0] = lattice[0];
		ss->lattice[1] = lattice[1];
		ss->lattice[2] = lattice[2];
		ss->lattice[3] = lattice[3];
		ss->lattice[4] = lattice[4];
		ss->lattice[5] = lattice[5];
		ss->lattice[6] = lattice[6];
		ss->lattice[7] = lattice[7];
		ss->lattice[8] = lattice[8];
	}
	return ss;
}

static void ShapeSet3_unfinalize(ShapeSet3 ss){
	if(ss->finalized){
		if(ss->use_bvh){
			BVH3_destroy(ss->bvh);
		}else{
			free(ss->v);
		}
		ss->finalized = 0;
	}
}

void ShapeSet3_destroy(ShapeSet3 ss){
	if(NULL == ss){ return; }
	ShapeSet3_unfinalize(ss);
	while(NULL != ss->lst){ // free the linked list
		shape3_list t = ss->lst;
		ss->lst = ss->lst->next;
		free(t);
	}
	free(ss);
}

int ShapeSet3_add(ShapeSet3 ss, const shape3 *s, int tag){
	unsigned int i, ilim = 1;
	static const int off[] = {
		0, 0, 0, // must be first
		 0,  0, -1,
		 1,  0, -1,
		-1,  0, -1,
		 0,  1, -1,
		 0, -1, -1,
		 1,  1, -1,
		 1, -1, -1,
		-1,  1, -1,
		-1, -1, -1,
		 1,  0, 0,
		-1,  0, 0,
		 0,  1, 0,
		 0, -1, 0,
		 1,  1, 0,
		 1, -1, 0,
		-1,  1, 0,
		-1, -1, 0,
		 0,  0, 1,
		 1,  0, 1,
		-1,  0, 1,
		 0,  1, 1,
		 0, -1, 1,
		 1,  1, 1,
		 1, -1, 1,
		-1,  1, 1,
		-1, -1, 1
	};
	if(NULL == ss){ return -1; }
	if(NULL == s){ return -2; }
	
	ShapeSet3_unfinalize(ss);
	
	if(ss->periodic){
		ilim += 26;
	}
	for(i = 0; i < ilim; ++i){
		shape3_list t = (shape3_list)malloc(sizeof(struct shape3_list_));
		memcpy(&(t->s), s, sizeof(shape3));
		SSDBG("Adding shape type=%d\n", t->s.type);
		t->s.o[0] += off[3*i+0] * ss->lattice[0] + off[3*i+1] * ss->lattice[3] + off[3*i+2] * ss->lattice[6];
		t->s.o[1] += off[3*i+0] * ss->lattice[1] + off[3*i+1] * ss->lattice[4] + off[3*i+2] * ss->lattice[7];
		t->s.o[2] += off[3*i+0] * ss->lattice[2] + off[3*i+1] * ss->lattice[5] + off[3*i+2] * ss->lattice[8];
		t->tag = tag;
		t->next = ss->lst;
		ss->lst = t;
		ss->n++;
	}
}


typedef struct{
	int index;
	const vecelem3 *v;
} shape3_iter_data;
static int shape3_iter(double b[6], int *tag, void *data){
	shape3_iter_data *d = (shape3_iter_data*)data;
	*tag = d->index; // the BVH's tag is actuall the index into the vector v
	b[0] = d->v[d->index].b.c[0] - d->v[d->index].b.h[0];
	b[1] = d->v[d->index].b.c[0] + d->v[d->index].b.h[0];
	b[2] = d->v[d->index].b.c[1] - d->v[d->index].b.h[1];
	b[3] = d->v[d->index].b.c[1] + d->v[d->index].b.h[1];
	b[4] = d->v[d->index].b.c[2] - d->v[d->index].b.h[2];
	b[5] = d->v[d->index].b.c[2] + d->v[d->index].b.h[2];
	d->index++;
	return 1;
}
static int vecelem3_cmp(const void *a, const void *b){
	const vecelem3 *e0 = (const vecelem3*)a;
	const vecelem3 *e1 = (const vecelem3*)b;
	return e1->tag - e0->tag;
}

static void ShapeSet3_finalize(ShapeSet3 ss){
	if(ss->finalized){ return; }
	
	ss->v = (vecelem3*)malloc(sizeof(vecelem3)*ss->n);
	// Insert linked list elements into vector v
	{
		shape3_list t = ss->lst;
		int index = 0;
		while(NULL != t){
			ss->v[index].s = &(t->s);
			shape3_get_aabb(&(t->s), &(ss->v[index].b));
			ss->v[index].tag = t->tag;
			++index;
			t = t->next;
		}
	}
	// Sort vector v so largest tags are at lower offsets
	qsort(ss->v, ss->n, sizeof(vecelem3), &vecelem3_cmp);

	if(ss->n >= ShapeSet3_threshold){
		shape3_iter_data d;
		d.index = 0;
		d.v = ss->v;
		ss->bvh = BVH3_new(ss->n, &shape3_iter, (void*)&d);
		ss->use_bvh = 1;
	}
	ss->finalized = 1;
}

static int query_pt3(int tag, const double b[6], void *data){
	int *best = (int*)data;
	if(tag > *best){ *best = tag; }
}

int ShapeSet3_query_pt(ShapeSet3 ss, const double p[3], const shape3 **s, int *tag){
	if(NULL == ss){ return -1; }
	if(NULL == p){ return -2; }
	
	ShapeSet3_finalize(ss);
	
	if(ss->use_bvh){
		int best = -1;
		BVH3_query_pt(ss->bvh, p, &query_pt3, &best);
		if(best >= 0){
			if(shape3_contains(ss->v[best].s, p)){
				if(NULL != tag){
					*tag = ss->v[best].tag;
				}
				if(NULL != s){
					*s = ss->v[best].s;
				}
				return 1;
			}
		}
	}else{
		int i;
		for(i = 0; i < ss->n; ++i){
			if(aabb3_contains(&(ss->v[i].b), p)){
				if(shape3_contains(ss->v[i].s, p)){
					if(NULL != tag){
						*tag = ss->v[i].tag;
					}
					if(NULL != s){
						*s = ss->v[i].s;
					}
					return 1;
				}
			}
		}
	}
	return 0;
}







/*

#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
#include <LuaCommon.h>

void ShapeSet_init(lua_State *L){
	static const struct luaL_reg ShapeSet2D_lib[] = {
		{"New", ShapeSet2D_New},
		{NULL,NULL}
	};
	static const struct luaL_Reg MaterialLibrary_obj[] = {
		{"Add", MaterialLibrary_Set},
		{NULL, NULL}
	};
}
*/