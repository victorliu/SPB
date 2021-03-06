#ifndef _LUAARG_H_INCLUDED_
#define _LUAARG_H_INCLUDED_

#include <lua.h>

/* recursive subtype */
#define luaarg_type_IGNORE 1
#define luaarg_type_TABLE 2
#define luaarg_type_LIST 3
/* scalar types */
#define luaarg_type_INT    4
#define luaarg_type_DOUBLE 5
#define luaarg_type_STRING 6 /* A char* pointer, use free() */
/* fixed size tuples */
#define luaarg_type_DOUBLE_VEC2   7
#define luaarg_type_DOUBLE_VEC3   8
#define luaarg_type_DOUBLE_MAT2X2 9
#define luaarg_type_DOUBLE_MAT3X3 10
#define luaarg_type_INT_VEC2   11
#define luaarg_type_INT_VEC3   12
#define luaarg_type_INT_MAT2X2 13
#define luaarg_type_INT_MAT3X3 14
#define luaarg_type_COMPLEX    15
#define luaarg_type_COMPLEX_CONSTITUTIVE_TENSOR3X3 16
/* variable size arrays */
#define luaarg_type_DOUBLE_VEC 17
#define luaarg_type_DOUBLE_MAT 18
#define luaarg_type_INT_VEC 19
#define luaarg_type_INT_MAT 20

typedef int (*luaarg_list_func)(lua_State *L, int keyind, int valind, void*);
typedef struct tag_luaarg_list_handler{
	luaarg_list_func func;
	void *data;
} luaarg_list_handler;

typedef struct tag_luaarg_double_vector{
	int n;
	double *v;
} luaarg_double_vector;

typedef struct tag_luaarg_double_matrix{
	int r,c;
	double *m;
} luaarg_double_matrix;

typedef struct tag_luaarg_int_vector{
	int n;
	int *v;
} luaarg_int_vector;

typedef struct tag_luaarg_int_matrix{
	int r,c;
	int *m;
} luaarg_int_matrix;

typedef struct tag_luaarg_argspec{
	const char *name;
	int type;
	int optional;
	void *valptr;
} luaarg_argspec;

void luaarg_parse(lua_State *L, int index, const luaarg_argspec *argspec);

#endif /* _LUAARG_H_INCLUDED_ */
